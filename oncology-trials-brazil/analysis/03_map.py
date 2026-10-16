"""Map of oncology trials per Brazilian state, plus a regional summary.

Reads data/processed/sites_brazil.csv (from 01_fetch_trials.py), assigns each site to a
state with clean_state() from 02_analyze.py and counts trials with at least one site in
each state. A trial with sites in several states counts once in each of them. State
boundaries come from the IBGE API (official 2022 state mesh) and are cached in data/raw/.

Run from the project folder:  python analysis/03_map.py
"""
import importlib.util
import json
import unicodedata
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Polygon

HERE = Path(__file__).resolve().parent
RAW = Path("data/raw")
SITES = Path("data/processed/sites_brazil.csv")
TRIALS = Path("data/processed/trials.csv")
RES = Path("results")
FIG = Path("figures")

MESH_URL = ("https://servicodados.ibge.gov.br/api/v3/malhas/paises/BR"
            "?formato=application/vnd.geo%2Bjson&intrarregiao=UF&qualidade=intermediaria")
STATES_URL = "https://servicodados.ibge.gov.br/api/v1/localidades/estados"

BINS = [0, 1, 10, 50, 200, 500, 1000, 2000]
BIN_LABELS = ["0", "1-9", "10-49", "50-199", "200-499", "500-999", "1,000+"]
COLOURS = ["#e6e6e6", "#d7ebf5", "#a6cde3", "#6baed6", "#3182bd", "#1b5d9b", "#08306b"]

# label offsets (lon, lat) for small states, so labels do not overlap
OFFSETS = {"DF": (2.6, 1.2), "RJ": (3.0, -2.0), "ES": (3.2, -0.6), "SE": (3.0, -1.2),
           "AL": (3.4, -0.4), "PE": (4.6, 0.3), "PB": (4.0, 1.0), "RN": (3.2, 1.6), "SC": (3.0, -0.6)}


def load_analyze_module():
    spec = importlib.util.spec_from_file_location("analyze", HERE / "02_analyze.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def plain(s):
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c)).lower()


def cached_json(url, path):
    if not path.exists():
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        path.write_bytes(r.content)
    return json.loads(path.read_text(encoding="utf-8"))


def rings(geometry):
    polys = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
    return [np.asarray(p[0]) for p in polys]   # outer rings only; states have no relevant holes


def ring_area_centroid(xy):
    x, y = xy[:, 0], xy[:, 1]
    cross = x[:-1] * y[1:] - x[1:] * y[:-1]
    a = cross.sum() / 2
    return abs(a), np.array([((x[:-1] + x[1:]) * cross).sum(), ((y[:-1] + y[1:]) * cross).sum()]) / (6 * a)


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    analyze = load_analyze_module()
    mesh = cached_json(MESH_URL, RAW / "ibge_states_mesh.geojson")
    states = pd.DataFrame([{"code": str(s["id"]), "uf": s["sigla"], "state": s["nome"],
                            "region": s["regiao"]["nome"]} for s in cached_json(STATES_URL, RAW / "ibge_states.json")])

    sites = pd.read_csv(SITES)
    n_trials = len(pd.read_csv(TRIALS))
    sites["state_clean"] = [analyze.clean_state(s, c) for s, c in zip(sites["state"], sites["city"])]
    name_to_uf = {plain(n): u for n, u in zip(states["state"], states["uf"])}
    sites["uf"] = sites["state_clean"].map(lambda s: name_to_uf.get(plain(s)))

    per_uf = sites.dropna(subset=["uf"]).drop_duplicates(["nct_id", "uf"]).groupby("uf").size()
    table = states.assign(trials=states["uf"].map(per_uf).fillna(0).astype(int))
    table["share_of_trials_pct"] = (100 * table["trials"] / n_trials).round(1)
    table = table.sort_values(["trials", "uf"], ascending=[False, True])
    table[["uf", "state", "region", "trials", "share_of_trials_pct"]].to_csv(RES / "trials_by_state_uf.csv", index=False)

    region = (sites.dropna(subset=["uf"]).merge(states[["uf", "region"]], on="uf")
              .drop_duplicates(["nct_id", "region"]).groupby("region").size()
              .reindex(["Sudeste", "Sul", "Nordeste", "Centro-Oeste", "Norte"]).fillna(0).astype(int))
    region_tbl = pd.DataFrame({"region": region.index, "trials": region.values,
                               "share_of_trials_pct": (100 * region.values / n_trials).round(1)})
    region_tbl.to_csv(RES / "trials_by_region.csv", index=False)
    n_partial = sites.loc[sites["uf"].isna(), "nct_id"].nunique()
    print(region_tbl.to_string(index=False))
    print(f"trials with at least one site whose state could not be identified: {n_partial}")

    cmap = ListedColormap(COLOURS)
    norm = BoundaryNorm(BINS, cmap.N)
    counts = table.set_index("code")["trials"]
    ufs = table.set_index("code")["uf"]

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"]})
    fig, (ax, axb) = plt.subplots(1, 2, figsize=(12, 7.2), width_ratios=[1.45, 1], layout="constrained")

    for feat in mesh["features"]:
        code = feat["properties"]["codarea"]
        n = int(counts.get(code, 0))
        best = (0, None)
        for ring in rings(feat["geometry"]):
            ax.add_patch(Polygon(ring, closed=True, facecolor=cmap(norm(n)), edgecolor="white", linewidth=0.6))
            area, c = ring_area_centroid(ring)
            best = max(best, (area, c), key=lambda t: t[0])
        uf = ufs[code]
        cx, cy = best[1]
        label = f"{uf}\n{n:,}"
        dark = n >= 200
        if uf in OFFSETS:
            dx, dy = OFFSETS[uf]
            ax.annotate(label, (cx, cy), (cx + dx, cy + dy), ha="center", va="center", fontsize=7.5,
                        arrowprops=dict(arrowstyle="-", color="#777777", lw=0.6))
        else:
            ax.text(cx, cy, label, ha="center", va="center", fontsize=7.5,
                    color="white" if dark else "#1a1a1a", fontweight="bold" if dark else "normal")
    ax.set_xlim(-74.5, -28.5)
    ax.set_ylim(-34.5, 6)
    ax.set_aspect(1 / np.cos(np.radians(15)))
    ax.axis("off")
    ax.set_title("Oncology trials with at least one site in the state", loc="left", fontsize=12, fontweight="bold")

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="#999999", lw=0.4) for c in COLOURS]
    ax.legend(handles, BIN_LABELS, title="Trials", loc="lower left", frameon=False, fontsize=8, title_fontsize=8.5)

    bars = table.iloc[::-1]
    axb.barh(bars["uf"], bars["trials"], color=[cmap(norm(v)) for v in bars["trials"]], edgecolor="#999999", lw=0.3)
    for y, v in enumerate(bars["trials"]):
        axb.text(v + 15, y, f"{v:,}", va="center", fontsize=7.5)
    axb.set_xlim(0, table["trials"].max() * 1.15)
    axb.set_xlabel("Trials with at least one site in the state")
    axb.tick_params(axis="y", labelsize=8)
    for side in ("top", "right"):
        axb.spines[side].set_visible(False)
    axb.set_title("Ranking", loc="left", fontsize=12, fontweight="bold")

    fig.text(0.01, -0.02,
             f"ClinicalTrials.gov snapshot of 2026-10-16: {n_trials:,} interventional oncology trials with a site in Brazil. "
             f"A trial with sites in several states is counted in each. {n_partial} trials have at least one site "
             "whose state could not be identified.\nStates assigned from the registry's state and city fields with "
             "the IBGE municipality list. Boundaries: IBGE state mesh.", fontsize=8, color="#444444")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"trials_map_brazil.{ext}", dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("map written to", FIG / "trials_map_brazil.png")


if __name__ == "__main__":
    main()
