"""Train and evaluate aqueous solubility (logS) models under two split strategies.

Splits (80/20, made once, before any preprocessing)
- random:   compounds assigned at random
- scaffold: Bemis-Murcko scaffold groups, largest first into train until 80% is
            reached; the rest (rarer chemotypes) form the test set. Every test
            molecule therefore has a ring system never seen in training.

Models
- mean:     predicts the training mean (floor any model must beat)
- esol:     linear model on the four ESOL inputs (Delaney 2004), refit on train
- rf:       random forest, RDKit descriptors + Morgan counts, fixed settings
- xgb:      gradient boosting, same features, depth / min_child_weight tuned by
            5-fold CV inside the training set (GroupKFold by scaffold for the
            scaffold split, so tuning sees the same kind of shift as the test)

Missing descriptor values are imputed with the training median inside each
pipeline, so test data never informs preprocessing. Test metrics carry 95%
bootstrap confidence intervals.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold, train_test_split
from sklearn.pipeline import make_pipeline
from xgboost import XGBRegressor

SEED = 42
N_BOOT = 2000
TEST_FRAC = 0.2

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

ESOL_COLS = ["esol_logp", "esol_mw", "esol_rotb", "esol_arom_prop"]
XGB_GRID = [{"max_depth": d, "min_child_weight": w} for d in (4, 6, 8) for w in (1, 5)]
XGB_FIXED = dict(n_estimators=800, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5,
                 tree_method="hist", random_state=SEED, n_jobs=-1)
SIM_BINS = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0001]


def rmse(y, p):
    return float(np.sqrt(mean_squared_error(y, p)))


def scaffold_split(scaffolds, test_frac):
    """Deterministic scaffold split: biggest scaffold groups go to train first."""
    groups = pd.Series(np.arange(len(scaffolds))).groupby(scaffolds.values).apply(list)
    groups = sorted(groups.values, key=lambda idx: (-len(idx), idx[0]))
    n_train_target = int(round((1 - test_frac) * len(scaffolds)))
    train, test = [], []
    for idx in groups:
        (train if len(train) + len(idx) <= n_train_target else test).extend(idx)
    return np.array(sorted(train)), np.array(sorted(test))


def bootstrap_ci(y, p, rng):
    n = len(y)
    stats = {"rmse": [], "r2": [], "mae": []}
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        stats["rmse"].append(rmse(y[i], p[i]))
        stats["r2"].append(r2_score(y[i], p[i]))
        stats["mae"].append(mean_absolute_error(y[i], p[i]))
    return {k: np.percentile(v, [2.5, 97.5]).round(3).tolist() for k, v in stats.items()}


def tune_xgb(X, y, groups):
    cv = GroupKFold(n_splits=5) if groups is not None else KFold(5, shuffle=True, random_state=SEED)
    rows = []
    for params in XGB_GRID:
        scores = []
        for tr, va in cv.split(X, y, groups):
            m = make_pipeline(SimpleImputer(strategy="median"), XGBRegressor(**XGB_FIXED, **params))
            m.fit(X[tr], y[tr])
            scores.append(rmse(y[va], m.predict(X[va])))
        rows.append({**params, "cv_rmse_mean": np.mean(scores), "cv_rmse_sd": np.std(scores)})
    cv_table = pd.DataFrame(rows).sort_values("cv_rmse_mean")
    best = cv_table.iloc[0][["max_depth", "min_child_weight"]].astype(int).to_dict()
    return best, cv_table


def nearest_train_similarity(smiles, train_idx, test_idx):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    train_fps = [fps[i] for i in train_idx]
    return np.array([max(DataStructs.BulkTanimotoSimilarity(fps[i], train_fps)) for i in test_idx])


def grouped_permutation_importance(model, X, y, feature_groups, rng, n_repeats=5):
    base = rmse(y, model.predict(X))
    rows = []
    for name, cols in feature_groups.items():
        deltas = []
        for _ in range(n_repeats):
            Xp = X.copy()
            Xp[:, cols] = Xp[rng.permutation(len(X))][:, cols]
            deltas.append(rmse(y, model.predict(Xp)) - base)
        rows.append({"feature": name, "rmse_increase_mean": np.mean(deltas),
                     "rmse_increase_sd": np.std(deltas)})
    return pd.DataFrame(rows).sort_values("rmse_increase_mean", ascending=False)


def main():
    rng = np.random.default_rng(SEED)
    meta = pd.read_csv(PROC / "molecules.csv", keep_default_na=False)
    meta["logS"] = meta.logS.astype(float)
    desc = pd.read_csv(PROC / "descriptors.csv")
    fp = np.load(PROC / "morgan_counts.npy").astype(np.float32)
    y = meta.logS.to_numpy()

    X_full = np.hstack([desc.to_numpy(np.float32), fp])
    X_esol = meta[ESOL_COLS].to_numpy(float)
    fp_cols = list(range(desc.shape[1], X_full.shape[1]))
    feature_groups = {c.removeprefix("rd_"): [i] for i, c in enumerate(desc.columns)}
    feature_groups["Morgan fingerprint (all 2048 bits)"] = fp_cols

    idx = np.arange(len(meta))
    tr_r, te_r = train_test_split(idx, test_size=TEST_FRAC, random_state=SEED)
    tr_s, te_s = scaffold_split(meta.scaffold, TEST_FRAC)
    splits = {"random": (tr_r, te_r), "scaffold": (tr_s, te_s)}

    metrics, preds, cv_tables, split_info = [], [], [], {}
    for split, (tr, te) in splits.items():
        scaf_tr, scaf_te = set(meta.scaffold.iloc[tr]), set(meta.scaffold.iloc[te])
        split_info[split] = {"n_train": len(tr), "n_test": len(te),
                             "test_scaffolds_seen_in_train": len(scaf_te & scaf_tr),
                             "test_scaffolds": len(scaf_te)}
        groups = meta.scaffold.iloc[tr].to_numpy() if split == "scaffold" else None

        best, cv_table = tune_xgb(X_full[tr], y[tr], groups)
        cv_tables.append(cv_table.assign(split=split))

        models = {
            "mean": None,
            "esol": make_pipeline(LinearRegression()),
            "rf": make_pipeline(SimpleImputer(strategy="median"),
                                RandomForestRegressor(n_estimators=500, max_features=0.3,
                                                      random_state=SEED, n_jobs=-1)),
            "xgb": make_pipeline(SimpleImputer(strategy="median"), XGBRegressor(**XGB_FIXED, **best)),
        }
        fitted = {}
        for name, model in models.items():
            Xtr, Xte = (X_esol[tr], X_esol[te]) if name == "esol" else (X_full[tr], X_full[te])
            if model is None:
                p = np.full(len(te), y[tr].mean())
            else:
                model.fit(Xtr, y[tr])
                p = model.predict(Xte)
                fitted[name] = model
            ci = bootstrap_ci(y[te], p, rng)
            metrics.append({"split": split, "model": name, "n_test": len(te),
                            "rmse": rmse(y[te], p), "rmse_ci": ci["rmse"],
                            "mae": mean_absolute_error(y[te], p), "mae_ci": ci["mae"],
                            "r2": r2_score(y[te], p), "r2_ci": ci["r2"],
                            "xgb_params": best if name == "xgb" else None})
            preds.append(pd.DataFrame({"split": split, "model": name, "row": te,
                                       "inchikey": meta.inchikey.iloc[te].to_numpy(),
                                       "observed": y[te], "predicted": p}))

        if split == "scaffold":
            sim = nearest_train_similarity(meta.smiles.tolist(), tr, te)
            p_xgb = fitted["xgb"].predict(X_full[te])
            b = pd.cut(sim, SIM_BINS, right=False)
            ad = pd.DataFrame({"sim": sim, "bin": b, "err": p_xgb - y[te]})
            ad_rows = []
            for interval, g in ad.groupby("bin", observed=True):
                e = g.err.to_numpy()
                boots = [np.sqrt(np.mean(rng.choice(e, len(e)) ** 2)) for _ in range(N_BOOT)]
                ad_rows.append({"similarity_bin": f"[{interval.left:.1f}, {min(interval.right, 1):.1f})",
                                "n": len(e), "rmse": float(np.sqrt(np.mean(e ** 2))),
                                "rmse_ci": np.percentile(boots, [2.5, 97.5]).round(3).tolist()})
            pd.DataFrame(ad_rows).to_csv(RESULTS / "error_by_similarity_scaffold_xgb.csv", index=False)
            pd.DataFrame({"row": te, "nn_tanimoto_to_train": sim}).to_csv(
                RESULTS / "test_similarity_scaffold.csv", index=False)

            imp = grouped_permutation_importance(fitted["xgb"], X_full[te], y[te], feature_groups, rng)
            imp.to_csv(RESULTS / "permutation_importance_scaffold_xgb.csv", index=False)

    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(RESULTS / "metrics.csv", index=False)
    pd.concat(preds).to_csv(RESULTS / "test_predictions.csv", index=False)
    pd.concat(cv_tables).to_csv(RESULTS / "xgb_cv_tuning.csv", index=False)
    (RESULTS / "split_info.json").write_text(json.dumps(split_info, indent=2))

    show = metrics_df[["split", "model", "n_test", "rmse", "rmse_ci", "mae", "r2", "r2_ci"]]
    print(show.round(3).to_string(index=False))
    print(json.dumps(split_info, indent=2))


if __name__ == "__main__":
    main()
