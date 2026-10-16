"""Curate AqSolDB structures and compute features.

Steps
1. Parse SMILES with RDKit and standardise: cleanup, keep the largest organic
   fragment (strips counter-ions), neutralise charges.
2. Keep organic molecules only: the parent has at least one carbon and only the
   elements in ALLOWED_ELEMENTS, and the original record holds no metal other than
   a common salt counter-ion. Metal complexes and inorganic salts are removed
   because the parent fragment does not describe their solubility.
3. Stripping salts can make two records the same parent molecule. Duplicates are
   merged by InChIKey: averaged if their logS values lie within MAX_DUP_RANGE,
   otherwise dropped as conflicting.
4. Features per molecule: the full RDKit 2D descriptor set, Morgan count
   fingerprints (radius 2, 2048 bits), the four ESOL inputs (Delaney 2004) and the
   Bemis-Murcko scaffold used for the scaffold split.

Every step only looks at one molecule at a time, so nothing here can leak
information between the later train and test sets.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "aqsoldb.csv"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

ALLOWED_ELEMENTS = {"C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se"}
# counter-ions that may be stripped from salts; any other metal in the original
# record (Cu, Zn, Cr, Co, ...) means a metal complex, which is excluded even when
# RDKit stores the metal as a separate fragment
SALT_COUNTERIONS = {"Na", "K", "Li", "Ca", "Mg"}
MAX_DUP_RANGE = 1.0  # log units
FP_BITS = 2048
FP_RADIUS = 2

_cleanup = rdMolStandardize.CleanupParameters()
_chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
_uncharger = rdMolStandardize.Uncharger()
_fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def standardise(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = rdMolStandardize.Cleanup(mol, _cleanup)
    mol = _chooser.choose(mol)
    mol = _uncharger.uncharge(mol)
    Chem.SanitizeMol(mol)
    return mol


def is_organic(parent, original):
    parent_elements = {a.GetSymbol() for a in parent.GetAtoms()}
    original_elements = {a.GetSymbol() for a in original.GetAtoms()}
    return ("C" in parent_elements and parent_elements <= ALLOWED_ELEMENTS
            and original_elements <= ALLOWED_ELEMENTS | SALT_COUNTERIONS)


def esol_inputs(mol):
    heavy = mol.GetNumHeavyAtoms()
    aromatic = sum(a.GetIsAromatic() for a in mol.GetAtoms())
    return {
        "esol_logp": Descriptors.MolLogP(mol),
        "esol_mw": Descriptors.MolWt(mol),
        "esol_rotb": Descriptors.NumRotatableBonds(mol),
        "esol_arom_prop": aromatic / heavy if heavy else 0.0,
    }


def main():
    PROC.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(RAW)
    log = {"raw_records": len(raw)}

    rows = []
    n_invalid = n_inorganic = 0
    for rec in raw.itertuples(index=False):
        mol = standardise(rec.SMILES)
        if mol is None or mol.GetNumHeavyAtoms() == 0:
            n_invalid += 1
            continue
        if not is_organic(mol, Chem.MolFromSmiles(rec.SMILES)):
            n_inorganic += 1
            continue
        rows.append({
            "source_id": rec.ID,
            "name": rec.Name,
            "smiles": Chem.MolToSmiles(mol),
            "inchikey": Chem.MolToInchiKey(mol),
            "logS": rec.Solubility,
            "source_sd": rec.SD,
            "source_occurrences": rec.Ocurrences,
            "changed_by_standardisation": Chem.MolToSmiles(mol) != Chem.CanonSmiles(rec.SMILES),
        })
    df = pd.DataFrame(rows)
    log.update({"invalid_smiles": n_invalid, "inorganic_or_metal": n_inorganic,
                "changed_by_standardisation": int(df.changed_by_standardisation.sum())})

    # merge parent-level duplicates
    grp = df.groupby("inchikey")
    spread = grp.logS.agg(lambda s: s.max() - s.min())
    sizes = grp.size()
    conflicting = spread[(sizes > 1) & (spread > MAX_DUP_RANGE)].index
    dup_keys = sizes[sizes > 1].index
    df = df[~df.inchikey.isin(conflicting)]
    agg = df.groupby("inchikey", as_index=False).agg(
        source_id=("source_id", "first"), name=("name", "first"), smiles=("smiles", "first"),
        logS=("logS", "mean"), source_sd=("source_sd", "max"),
        source_occurrences=("source_occurrences", "sum"), changed_by_standardisation=("changed_by_standardisation", "any"),
        n_merged=("logS", "size"),
    )
    log.update({"parent_duplicate_groups": len(dup_keys),
                "conflicting_groups_dropped": len(conflicting),
                "final_molecules": len(agg)})

    # features
    desc_rows, fp_rows, extra = [], [], []
    for smi in agg.smiles:
        mol = Chem.MolFromSmiles(smi)
        desc_rows.append(Descriptors.CalcMolDescriptors(mol))
        fp_rows.append(_fpgen.GetCountFingerprintAsNumPy(mol).astype(np.uint16))
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        extra.append({**esol_inputs(mol), "scaffold": scaffold})

    desc = pd.DataFrame(desc_rows).add_prefix("rd_").astype(float)
    # Ipc grows exponentially with molecular size (up to ~1e54 here), beyond float32
    desc["rd_Ipc"] = np.log10(desc["rd_Ipc"].clip(lower=1.0))
    desc = desc.rename(columns={"rd_Ipc": "rd_log10_Ipc"})
    desc = desc.replace([np.inf, -np.inf], np.nan)
    meta = pd.concat([agg.reset_index(drop=True), pd.DataFrame(extra)], axis=1)

    meta.to_csv(PROC / "molecules.csv", index=False)
    desc.to_csv(PROC / "descriptors.csv", index=False)
    np.save(PROC / "morgan_counts.npy", np.vstack(fp_rows))

    log.update({
        "rdkit_descriptors": desc.shape[1],
        "descriptor_cells_missing": int(desc.isna().sum().sum()),
        "acyclic_molecules": int((meta.scaffold == "").sum()),
        "unique_scaffolds": int(meta.scaffold.nunique()),
    })
    multi = raw[raw.Ocurrences > 1].SD
    log["noise_floor"] = {
        "compounds_with_multiple_sources": int(len(multi)),
        "median_between_source_sd": round(float(multi.median()), 3),
        "rms_between_source_sd": round(float(np.sqrt((multi ** 2).mean())), 3),
    }
    (RESULTS / "curation_log.json").write_text(json.dumps(log, indent=2))
    print(json.dumps(log, indent=2))


if __name__ == "__main__":
    main()
