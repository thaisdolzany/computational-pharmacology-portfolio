"""Download AqSolDB (Sorkun et al. 2019, Sci Data 6:143) from Harvard Dataverse.

The file is checked against the MD5 published by Dataverse, and the download date
is written next to it so the run can be traced later.
"""
import hashlib
import json
from datetime import date
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

# Dataverse dataset doi:10.7910/DVN/OVHAW8, version 1 (CC0 1.0)
FILE_ID = 3407241
EXPECTED_MD5 = "02f732e52de921b8adaf4ef03cf747bf"
URL = f"https://dataverse.harvard.edu/api/access/datafile/{FILE_ID}?format=original"
OUT = RAW / "aqsoldb.csv"


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    if not OUT.exists():
        resp = requests.get(URL, timeout=120)
        resp.raise_for_status()
        OUT.write_bytes(resp.content)

    md5 = hashlib.md5(OUT.read_bytes()).hexdigest()
    if md5 != EXPECTED_MD5:
        raise RuntimeError(f"MD5 mismatch for {OUT.name}: got {md5}, expected {EXPECTED_MD5}")

    log = {
        "source": "AqSolDB, Harvard Dataverse doi:10.7910/DVN/OVHAW8 (v1, CC0 1.0)",
        "file_id": FILE_ID,
        "md5": md5,
        "downloaded": date.today().isoformat(),
    }
    (RAW / "download_log.json").write_text(json.dumps(log, indent=2))
    print(f"saved {OUT.name} ({OUT.stat().st_size / 1e6:.1f} MB), md5 ok")


if __name__ == "__main__":
    main()
