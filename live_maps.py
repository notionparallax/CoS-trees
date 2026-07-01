"""
Compare 2017 FOI data with the live City of Sydney open data portal.

Run with:  .venv/Scripts/python live_maps.py
Outputs:   out/ directory — PNG maps using current data
"""

import os
from random import seed, shuffle

import geodatasets
import geopandas
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point, box

from tree_common_to_latin_names_map import tree_common_names

os.makedirs("out", exist_ok=True)

cm_to_inch = 2.54
plt.rcParams["figure.figsize"] = [39 / cm_to_inch, 49 / cm_to_inch]

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print("Loading live CoS data from ArcGIS open data portal…")
url = "https://opendata.arcgis.com/datasets/cityofsydney::trees.geojson"
gdf = geopandas.read_file(url)

print(f"\n{len(gdf):,} trees loaded.")
print(f"CRS: {gdf.crs}")
print(f"\nColumns:\n{list(gdf.columns)}\n")
print(gdf.dtypes)
print(f"\nSample:\n{gdf.sample(3).to_string()}")

# Normalise column names to lowercase for consistency
gdf.columns = [c.lower() for c in gdf.columns]

# Ensure a 'species' column exists.
# After lowercasing, ArcGIS dataset uses 'speciesname'.
species_candidates = ["species", "speciesname", "botanical_name", "botanicalname", "scientific_name"]
species_col = next((c for c in species_candidates if c in gdf.columns), None)
if species_col and species_col != "species":
    gdf = gdf.rename(columns={species_col: "species"})
    print(f"\nRenamed '{species_col}' -> 'species'")
elif species_col is None:
    print("\nWARNING: no species column found; available columns printed above.")
    gdf["species"] = "unknown"

gdf["species"] = [x if isinstance(x, str) else "unknown" for x in gdf["species"]]

# ArcGIS GeoJSON endpoints often serve data in the layer's native projected CRS
# rather than the GeoJSON-standard WGS84. Detect this by checking coordinate range.
sample = gdf.geometry.iloc[0]
if abs(sample.x) > 1000 or abs(sample.y) > 1000:
    # Coordinates look projected (MGA56/EPSG:28356 for Sydney)
    if gdf.crs is None or gdf.crs.to_epsg() != 28356:
        gdf = gdf.set_crs(epsg=28356, allow_override=True)
        print("\nProjected coordinates detected — assigned EPSG:28356 (MGA56)")
    gdf = gdf.to_crs(epsg=4326)
elif gdf.crs is not None and gdf.crs.to_epsg() != 4326:
    gdf = gdf.to_crs(epsg=4326)

print(f"CRS after normalisation: {gdf.crs}")

# ---------------------------------------------------------------------------
# Base geography — suburb outlines clipped to tree extent
# ---------------------------------------------------------------------------

# Compute tight bounds from the actual tree data with a small padding
PAD = 0.005  # ~500 m
minx, miny, maxx, maxy = gdf.total_bounds
xlim = (minx - PAD, maxx + PAD)
ylim = (miny - PAD, maxy + PAD)

aus_poas = geopandas.read_file("geopandas-blog/aus_poas.shp")
bbox_geom = box(xlim[0], ylim[0], xlim[1], ylim[1])
syd = aus_poas[aus_poas.intersects(bbox_geom)]

# ---------------------------------------------------------------------------
# Map 1 — all trees
# ---------------------------------------------------------------------------

print("\nRendering map: all trees…")
fig, ax = plt.subplots()
syd.plot(ax=ax, color="ghostwhite", edgecolor="black")
gdf.plot(ax=ax, color="darkgreen", marker=".", markersize=1, alpha=0.15)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
ax.set_title(f"All trees in the City of Sydney — live open data ({len(gdf):,} trees)")
plt.tight_layout()
plt.savefig("out/all_trees_live.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved out/all_trees_live.png")

# ---------------------------------------------------------------------------
# Map 2 — top species bar chart
# ---------------------------------------------------------------------------

print("Rendering chart: top species…")
top = 30
fig, ax = plt.subplots(figsize=(14, 6))
gdf["species"].value_counts()[:top].plot(kind="bar", ax=ax, rot=90)
ax.set_title(f"Top {top} most common trees — City of Sydney (live data)")
ax.set_ylabel("Count")
plt.tight_layout()
plt.savefig("out/top_species_live.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved out/top_species_live.png")

# ---------------------------------------------------------------------------
# Map 3 — figs
# ---------------------------------------------------------------------------

print("Rendering map: figs…")
figs = {
    "Ficus microcarpa var hillii": {"colour": "C0", "marker": "+", "alpha": 0.3},
    "Ficus benjamina":             {"colour": "C1", "marker": ".", "alpha": 0.3},
    "Ficus rubiginosa":            {"colour": "C2", "marker": "1", "alpha": 0.3},
    "Ficus macrophylla":           {"colour": "C3", "marker": "*", "alpha": 0.9},
    "Ficus sp.":                   {"colour": "C4", "marker": "p", "alpha": 0.3},
    "Ficus carica":                {"colour": "C5", "marker": "o", "alpha": 0.3},
    "Ficus elastica":              {"colour": "C6", "marker": "s", "alpha": 0.3},
}
fig_gdf = gdf[[("Ficus" in x) for x in gdf["species"]]]
print(f"  {len(fig_gdf)} figs found")

fig_map, ax = plt.subplots()
syd.plot(ax=ax, color="ghostwhite", edgecolor="black")
for species, style in figs.items():
    sub = fig_gdf[fig_gdf["species"] == species]
    if len(sub):
        sub.plot(ax=ax, color=style["colour"], marker=style["marker"],
                 alpha=style["alpha"], markersize=4)
ax.legend(
    handles=[mpatches.Patch(color=figs[f]["colour"], label=f) for f in figs],
    loc="upper left", fontsize=7,
)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
ax.set_title("Figs in the City of Sydney — live open data")
plt.tight_layout()
plt.savefig("out/figs_live.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved out/figs_live.png")

# ---------------------------------------------------------------------------
# Map 4 — all species by colour/marker (equiv. to the original all_tree_map)
# ---------------------------------------------------------------------------

print("Rendering map: all species…")
markers = [".", "o", "v", "^", "<", ">", "1", "2", "3", "4",
           "8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|", "_"]
tab = list(mcolors.TABLEAU_COLORS.keys())
base = list(mcolors.BASE_COLORS.keys())
tab.extend(base)
colours = [c for c in tab if c != "w"]
idents = [{"m": m, "c": c} for m in markers for c in colours]
seed(42)
shuffle(idents)

tree_counts = gdf["species"].value_counts()

fig_all, ax = plt.subplots()
syd.plot(ax=ax, color="ghostwhite", edgecolor="white")
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)

legend_handles = []
for i, (tree, count) in enumerate(tree_counts.items()):
    ident = idents[i % len(idents)]
    sub = gdf[gdf["species"] == tree]
    try:
        common = "—" + tree_common_names[tree][0].title()
    except (KeyError, IndexError):
        common = ""
    label = rf"$\it{{{tree}}}${common} ({count})" if count >= 50 else None
    sub.plot(ax=ax, color=ident["c"], marker=ident["m"],
             alpha=0.4, markersize=3, label=label)
    if label:
        legend_handles.append(
            mpatches.Patch(color=ident["c"],
                           label=rf"$\it{{{tree}}}${common} ({count})")
        )

ax.legend(handles=legend_handles, loc="upper right",
          markerscale=2, prop={"size": 5}, labelspacing=0)
ax.tick_params(labelsize=6)
title = f"{len(gdf):,} Trees in the City of Sydney (live data)"
ax.set_title(title)
plt.tight_layout()
plt.savefig("out/all_tree_map_live.png", dpi=300, bbox_inches="tight")
plt.close()
print("  saved out/all_tree_map_live.png")

# ---------------------------------------------------------------------------
# Map 5 — new: tree height (only possible with live data)
# ---------------------------------------------------------------------------

height_col = next((c for c in ["treeheight", "height", "tree_height", "height_m"] if c in gdf.columns), None)
if height_col:
    print(f"Rendering map: tree height (column '{height_col}')…")
    height_gdf = gdf[gdf[height_col].notna()].copy()
    height_gdf[height_col] = pd.to_numeric(height_gdf[height_col], errors="coerce")
    height_gdf = height_gdf[height_gdf[height_col] > 0]

    fig_h, ax = plt.subplots()
    syd.plot(ax=ax, color="ghostwhite", edgecolor="black")
    # Cap the colour scale at the 95th percentile so outliers don't wash out the scale
    vmax = height_gdf[height_col].quantile(0.95)
    height_gdf.plot(ax=ax, column=height_col, cmap="YlGn",
                    markersize=2, alpha=0.6, legend=True, vmin=0, vmax=vmax,
                    legend_kwds={"label": f"Height (m, capped at {vmax:.0f}m / 95th pctile)", "shrink": 0.6})
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title("Tree height — City of Sydney (live data)")
    plt.tight_layout()
    plt.savefig("out/tree_height_live.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved out/tree_height_live.png")
else:
    print("  (no height column found — skipping height map)")

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

print("\nAll done.  Outputs in out/")
print("\nLive dataset column summary:")
print(gdf.describe(include="all").T[["count", "unique", "top"]].to_string())
