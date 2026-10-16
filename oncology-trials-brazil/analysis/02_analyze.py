"""Summarise the Brazilian oncology trial landscape from the files written by
01_fetch_trials.py. Writes tables to results/ and figures to figures/.

Run from the project folder:  python analysis/02_analyze.py
"""

import re
import unicodedata
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

TRIALS = Path("data/processed/trials.csv")
SITES = Path("data/processed/sites_brazil.csv")
RES = Path("results")
FIG = Path("figures")

# First match wins, so the more specific terms come first.
TUMOUR_TERMS = [
    ("Breast", r"breast"),
    ("Lung", r"lung|nsclc|sclc"),
    ("Prostate", r"prostat"),
    ("Colorectal", r"colorectal|colon|rectal"),
    ("Melanoma", r"melanoma"),
    ("Lymphoma", r"lymphoma"),
    ("Leukaemia", r"leuk"),
    ("Myeloma", r"myeloma"),
    ("Gastric/oesophageal", r"gastric|stomach|esophag|oesophag"),
    ("Liver", r"hepatocellular|liver"),
    ("Kidney/bladder", r"renal|kidney|urothelial|bladder"),
    ("Ovarian/cervical/endometrial", r"ovar|cervi|endometri|uter"),
    ("Head and neck", r"head and neck|oropharyn|laryn|nasopharyn"),
    ("Pancreas", r"pancrea"),
    ("Brain", r"glioma|glioblastoma|brain"),
    ("Sarcoma", r"sarcoma"),
]

STATE_NAMES = {
    "AC": "Acre", "AL": "Alagoas", "AP": "Amapa", "AM": "Amazonas", "BA": "Bahia",
    "CE": "Ceara", "DF": "Distrito Federal", "ES": "Espirito Santo", "GO": "Goias",
    "MA": "Maranhao", "MT": "Mato Grosso", "MS": "Mato Grosso do Sul",
    "MG": "Minas Gerais", "PA": "Para", "PB": "Paraiba", "PR": "Parana",
    "PE": "Pernambuco", "PI": "Piaui", "RJ": "Rio de Janeiro",
    "RN": "Rio Grande do Norte", "RS": "Rio Grande do Sul", "RO": "Rondonia",
    "RR": "Roraima", "SC": "Santa Catarina", "SP": "Sao Paulo", "SE": "Sergipe",
    "TO": "Tocantins",
}
NAME_TO_STATE = {v.lower(): v for v in STATE_NAMES.values()}


def plain(text):
    text = str(text)
    try:
        text = text.encode("latin-1").decode("utf-8")  # undo mojibake such as "SÃ£o"
    except (UnicodeEncodeError, UnicodeDecodeError):
        pass
    text = unicodedata.normalize("NFKD", text)
    return "".join(c for c in text if not unicodedata.combining(c)).strip()


# Registry entries spell states in many ways (typos, abbreviations, city
# names in the state field). Each pattern is matched against the lower-case,
# accent-free text.
STATE_PATTERNS = [
    ("Sao Paulo", r"s\w*o paulo|\bsp\b|santo andre|paulo prata|de ao paulo"),
    ("Rio Grande do Sul", r"rio gran\w* d\w* s|\brs\b|porto? alegre|porte alegre"),
    ("Rio de Janeiro", r"rio d\w* jan|\brj\b"),
    ("Parana", r"parana|curitiba|\bpr\b"),
    ("Minas Gerais", r"minas|\bmg\b|\bbh\b"),
    ("Santa Catarina", r"santa cat"),
    ("Bahia", r"bahia"),
    ("Ceara", r"ceara|fortal"),
    ("Goias", r"goias"),
    ("Distrito Federal", r"federal district|distrito federal|brasilia"),
    ("Espirito Santo", r"espirito|vitoria"),
    ("Rondonia", r"rond"),
]

# Used when the state field is empty: the cities that host most sites.
CITY_TO_STATE = {
    "sao paulo": "Sao Paulo", "barretos": "Sao Paulo", "santo andre": "Sao Paulo",
    "sao jose do rio preto": "Sao Paulo", "ribeirao preto": "Sao Paulo",
    "jau": "Sao Paulo", "campinas": "Sao Paulo", "sorocaba": "Sao Paulo",
    "sao caetano do sul": "Sao Paulo", "taubate": "Sao Paulo", "botucatu": "Sao Paulo",
    "porto alegre": "Rio Grande do Sul", "ijui": "Rio Grande do Sul",
    "caxias do sul": "Rio Grande do Sul", "passo fundo": "Rio Grande do Sul",
    "santa maria": "Rio Grande do Sul", "pelotas": "Rio Grande do Sul",
    "lajeado": "Rio Grande do Sul", "rio grande": "Rio Grande do Sul",
    "santa cruz do sul": "Rio Grande do Sul",
    "rio de janeiro": "Rio de Janeiro", "curitiba": "Parana", "londrina": "Parana",
    "belo horizonte": "Minas Gerais", "ipatinga": "Minas Gerais", "uberlandia": "Minas Gerais",
    "salvador": "Bahia", "natal": "Rio Grande do Norte", "fortaleza": "Ceara",
    "florianopolis": "Santa Catarina", "itajai": "Santa Catarina",
    "joinville": "Santa Catarina", "blumenau": "Santa Catarina",
    "goiania": "Goias", "brasilia": "Distrito Federal", "vitoria": "Espirito Santo",
    "cachoeiro de itapemirim": "Espirito Santo", "recife": "Pernambuco",
    "teresina": "Piaui", "belem": "Para",
}


def clean_state(state, city=None):
    if not pd.isna(state):
        v = plain(state)
        if v.upper() in STATE_NAMES:
            return STATE_NAMES[v.upper()]
        if v.lower() in NAME_TO_STATE:
            return NAME_TO_STATE[v.lower()]
        low = re.sub(r"^(state of|estado d[eo])\s+", "", v.lower())
        if low in NAME_TO_STATE:
            return NAME_TO_STATE[low]
        for name, pattern in STATE_PATTERNS:
            if re.search(pattern, low):
                return name
    if city is not None and not pd.isna(city):
        found = CITY_TO_STATE.get(plain(city).lower())
        if found:
            return found
    return "Not reported"


def tumour_type(conditions):
    text = str(conditions).lower()
    for label, pattern in TUMOUR_TERMS:
        if re.search(pattern, text):
            return label
    return "Other or mixed"


def phase_label(phases):
    if pd.isna(phases) or phases in ("", "NA"):
        return "Not applicable"
    parts = [p.replace("EARLY_PHASE1", "Early 1").replace("PHASE", "") for p in phases.split(";")]
    return "Phase " + "/".join(parts)


def counts(series):
    """value_counts with ties broken by name, so the order does not depend on the pandas version."""
    vc = series.value_counts()
    return vc.rename_axis("label").reset_index(name="count").sort_values(
        ["count", "label"], ascending=[False, True]).set_index("label")["count"].rename_axis(series.name)


def save_bar(series, title, xlabel, path, top=None):
    s = series.head(top) if top else series
    fig, ax = plt.subplots(figsize=(7, max(3, 0.35 * len(s))))
    s.iloc[::-1].plot.barh(ax=ax, color="#0072B2")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    RES.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    trials = pd.read_csv(TRIALS)
    sites = pd.read_csv(SITES)

    trials["start_year"] = pd.to_datetime(trials["start_date"], errors="coerce").dt.year
    trials["sponsor_type"] = trials["sponsor_class"].map({"INDUSTRY": "Industry"}).fillna("Academic or public")
    trials["phase"] = trials["phases"].apply(phase_label)
    trials["tumour"] = trials["conditions"].apply(tumour_type)
    trials["multinational"] = trials["n_countries"] > 1

    summary = {
        "trials": len(trials),
        "industry_share": round((trials["sponsor_type"] == "Industry").mean(), 3),
        "multinational_share": round(trials["multinational"].mean(), 3),
        "recruiting_now": int((trials["status"] == "RECRUITING").sum()),
        "brazilian_sites": len(sites),
        "median_brazilian_sites_per_trial": float(trials["n_sites_brazil"].median()),
    }
    pd.Series(summary).to_csv(RES / "summary.csv", header=["value"])

    by_year = (
        trials.dropna(subset=["start_year"])
        .query("2010 <= start_year <= 2026")
        .groupby(["start_year", "sponsor_type"]).size().unstack(fill_value=0)
    )
    by_year.to_csv(RES / "trials_by_start_year.csv")
    fig, ax = plt.subplots(figsize=(8, 4))
    styles = {"Industry": ("#D55E00", "o"), "Academic or public": ("#0072B2", "s")}
    for col in by_year.columns:
        colour, marker = styles.get(col, ("#555555", "^"))
        ax.plot(by_year.index, by_year[col], marker=marker, color=colour, label=col)
    years = by_year.index.astype(int)
    last = years.max()
    top = by_year.to_numpy().max() * 1.15
    ax.set_ylim(0, top)
    ax.axvspan(last - 0.5, last + 0.5, color="#eeeeee", zorder=0)
    ax.text(last, top * 0.97, "partial\nyear", ha="center", va="top", fontsize=8, color="#555555")
    ax.set_xticks(range(years.min(), last + 1, 2))
    ax.legend(title="Lead sponsor", frameon=False, loc="upper left")
    ax.set_title("Oncology trials with a site in Brazil, by start year")
    ax.set_xlabel("Start year")
    ax.set_ylabel("Trials started")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "trials_by_start_year.png", dpi=200)
    plt.close(fig)

    phase = counts(trials["phase"])
    phase.to_csv(RES / "trials_by_phase.csv")
    save_bar(phase, "Trials by phase", "Trials", FIG / "trials_by_phase.png")

    tumour = counts(trials["tumour"])
    tumour.to_csv(RES / "trials_by_tumour.csv")
    save_bar(tumour, "Trials by tumour type", "Trials", FIG / "trials_by_tumour.png")

    sponsors = counts(trials.loc[trials["sponsor_type"] == "Industry", "lead_sponsor"])
    sponsors.to_csv(RES / "industry_sponsors.csv")
    save_bar(sponsors, "Top 15 industry sponsors", "Trials", FIG / "industry_sponsors.png", top=15)

    sites["state_clean"] = [clean_state(s, c) for s, c in zip(sites["state"], sites["city"])]
    states = counts(sites.drop_duplicates(["nct_id", "state_clean"])["state_clean"])
    states.to_csv(RES / "trials_by_state.csv")
    save_bar(states, "Trials with at least one site, by state", "Trials", FIG / "trials_by_state.png", top=15)

    cities = counts(sites.drop_duplicates(["nct_id", "city"])["city"])
    cities.to_csv(RES / "trials_by_city.csv")

    status = counts(trials["status"])
    status.to_csv(RES / "trials_by_status.csv")

    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
