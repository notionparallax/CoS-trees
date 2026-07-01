"""
Fetch the USYD CampusFlora plant location data and save it to a local CSV.

Run this whenever you want to refresh the cached data:
    .venv/Scripts/python update_usyd_data.py

The main live_maps.py script reads from the CSV so it doesn't hit their
server on every run.
"""

import requests
import pandas as pd
from pathlib import Path

URL = "https://campusflora.sydney.edu.au/species.json"
OUT = Path("in/usyd_campus_flora.csv")

print(f"Fetching {URL} …")
resp = requests.get(URL, timeout=30)
resp.raise_for_status()

rows = []
for sp in resp.json():
    for loc in sp.get("species_locations", []):
        if not loc.get("removed", True):
            rows.append({
                "species": sp["genus_species"],
                "longitude": float(loc["longitude"]),
                "latitude": float(loc["latitude"]),
                "source": "University of Sydney",
            })

df = pd.DataFrame(rows)
OUT.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT, index=False)
print(f"Saved {len(df):,} plant locations ({df['species'].nunique()} species) → {OUT}")
