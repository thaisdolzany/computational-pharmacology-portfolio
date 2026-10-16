"""Download interventional cancer trials with at least one site in Brazil
from the ClinicalTrials.gov API v2 and save one row per trial plus one row
per Brazilian site.

Run from the project folder:  python analysis/01_fetch_trials.py
"""

import json
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://clinicaltrials.gov/api/v2/studies"
OUT_RAW = Path("data/raw/trials_raw.jsonl")
OUT_TRIALS = Path("data/processed/trials.csv")
OUT_SITES = Path("data/processed/sites_brazil.csv")

PARAMS = {
    "query.cond": "cancer OR neoplasm OR tumor OR carcinoma OR lymphoma OR leukemia OR melanoma OR myeloma OR sarcoma",
    "query.locn": "Brazil",
    "filter.advanced": "AREA[StudyType]INTERVENTIONAL",
    "pageSize": 1000,
    "format": "json",
}


def fetch_all():
    studies, token, page = [], None, 0
    while True:
        params = dict(PARAMS)
        if token:
            params["pageToken"] = token
        r = requests.get(API, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        batch = data.get("studies", [])
        studies.extend(batch)
        page += 1
        print(f"page {page}: {len(batch)} studies (total {len(studies)})")
        token = data.get("nextPageToken")
        if not token:
            return studies
        time.sleep(0.5)


def get(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


def flatten(study):
    p = study.get("protocolSection", {})
    locations = get(p, "contactsLocationsModule", "locations", default=[])
    countries = {loc.get("country") for loc in locations if loc.get("country")}
    trial = {
        "nct_id": get(p, "identificationModule", "nctId"),
        "title": get(p, "identificationModule", "briefTitle"),
        "status": get(p, "statusModule", "overallStatus"),
        "start_date": get(p, "statusModule", "startDateStruct", "date"),
        "phases": ";".join(get(p, "designModule", "phases", default=[])),
        "enrollment": get(p, "designModule", "enrollmentInfo", "count"),
        "lead_sponsor": get(p, "sponsorCollaboratorsModule", "leadSponsor", "name"),
        "sponsor_class": get(p, "sponsorCollaboratorsModule", "leadSponsor", "class"),
        "conditions": ";".join(get(p, "conditionsModule", "conditions", default=[])),
        "n_countries": len(countries),
        "n_sites_total": len(locations),
        "n_sites_brazil": sum(1 for loc in locations if loc.get("country") == "Brazil"),
    }
    sites = [
        {
            "nct_id": trial["nct_id"],
            "facility": loc.get("facility"),
            "city": loc.get("city"),
            "state": loc.get("state"),
        }
        for loc in locations
        if loc.get("country") == "Brazil"
    ]
    return trial, sites


def main():
    for path in (OUT_RAW, OUT_TRIALS, OUT_SITES):
        path.parent.mkdir(parents=True, exist_ok=True)
    studies = fetch_all()
    with OUT_RAW.open("w", encoding="utf-8") as fh:
        for s in studies:
            fh.write(json.dumps(s) + "\n")
    trials, sites = [], []
    for s in studies:
        t, st = flatten(s)
        trials.append(t)
        sites.extend(st)
    trials = pd.DataFrame(trials)
    trials = trials[trials["n_sites_brazil"] > 0]
    pd.DataFrame(sites).to_csv(OUT_SITES, index=False)
    trials.to_csv(OUT_TRIALS, index=False)
    print(f"saved {len(trials)} trials and {len(sites)} Brazilian sites")


if __name__ == "__main__":
    main()
