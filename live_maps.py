"""
Compare 2017 FOI data with the live City of Sydney open data portal.

Run with:  .venv/Scripts/python live_maps.py
Outputs:   out/ directory — PNG maps using current data
"""

import os
from random import seed, shuffle

import requests

import geodatasets
import geopandas
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point, box

from tree_common_to_latin_names_map import tree_common_names

os.makedirs("out", exist_ok=True)


def clean_ax(ax):
    """Remove spines, ticks and tick labels from a map axes."""
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


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
gdf["source"] = "City of Sydney (live)"

# ---------------------------------------------------------------------------
# Load supplementary data: Ryde (2013) and RBG / Centennial Parklands
# ---------------------------------------------------------------------------

print("\nLoading Ryde tree data (2013)…")
ryde_gdf = geopandas.read_file(os.path.join("in", "ryde", "public-trees-2013.shp"))
ryde_gdf = ryde_gdf.set_crs(epsg=28356).to_crs(epsg=4326)
ryde_gdf["species"] = "unknown"
ryde_gdf["source"] = "City of Ryde (2013)"
print(f"  {len(ryde_gdf):,} trees")

print("Loading RBG / Centennial Parklands tree data…")
rbg_frames = []
for file_name, short_name in [
    ("BGCP-AustralianBotanicGardens-MtAnnanTreeList.xlsx", "Mt Annan"),
    ("BGCP-BlueMountainsBotanicGardens-TomahTreeList.xlsx", "Mt Tomah"),
    ("BGCP-CentennialParklandsTreeList.xlsx",               "Centennial Parklands"),
    ("BGCP-RoyalBotanicGarden-SydneyTreeList.xlsx",         "Royal Botanic Garden"),
]:
    path = os.path.join("in", "RBG", file_name)
    df = pd.read_excel(path, skiprows=3)
    df = df[df["ItemCoordLongDD"].notna() & df["ItemCoordLatDD"].notna()]
    df["species"] = df["GenusSpecies"]
    df["source"] = short_name
    rbg_frames.append(df)
rbg_df = pd.concat(rbg_frames, ignore_index=True)
rbg_gdf = geopandas.GeoDataFrame(
    rbg_df,
    geometry=[Point(xy) for xy in zip(rbg_df["ItemCoordLongDD"], rbg_df["ItemCoordLatDD"])],
    crs="EPSG:4326",
)
print(f"  {len(rbg_gdf):,} trees ({', '.join(rbg_df['source'].unique())})")

print("Loading USYD CampusFlora data…")
_usyd_resp = requests.get("https://campusflora.sydney.edu.au/species.json", timeout=30)
_usyd_resp.raise_for_status()
_usyd_rows = []
for _sp in _usyd_resp.json():
    for _loc in _sp.get("species_locations", []):
        if not _loc.get("removed", True):
            _usyd_rows.append({
                "species": _sp["genus_species"],
                "longitude": float(_loc["longitude"]),
                "latitude": float(_loc["latitude"]),
                "source": "University of Sydney",
            })
_usyd_df = pd.DataFrame(_usyd_rows)
usyd_gdf = geopandas.GeoDataFrame(
    _usyd_df,
    geometry=geopandas.points_from_xy(_usyd_df["longitude"], _usyd_df["latitude"]),
    crs="EPSG:4326",
)
print(f"  {len(usyd_gdf):,} plant locations ({_usyd_df['species'].nunique()} species)")

# ---------------------------------------------------------------------------
# Filter to within 10 km of Sydney CBD; build combined GDF
# ---------------------------------------------------------------------------

CBD_LON, CBD_LAT = 151.207, -33.869
cbd_proj = (
    geopandas.GeoDataFrame(geometry=[Point(CBD_LON, CBD_LAT)], crs="EPSG:4326")
    .to_crs(epsg=28356)
    .geometry.iloc[0]
)
RADIUS_TIGHT = 4_000   # 4 km  — covers CoS, RBG, Centennial; excludes Ryde
RADIUS_WIDE  = 20_000  # 20 km — covers all of Ryde and beyond

def within_radius(df, radius):
    return df[df.to_crs(epsg=28356).geometry.distance(cbd_proj) <= radius]

ryde_near = within_radius(ryde_gdf, RADIUS_TIGHT)
rbg_near  = within_radius(rbg_gdf,  RADIUS_TIGHT)
usyd_near = within_radius(usyd_gdf, RADIUS_TIGHT)
ryde_wide = within_radius(ryde_gdf, RADIUS_WIDE)

print(f"\nWithin {RADIUS_TIGHT/1000:.0f} km of CBD (tight):")
print(f"  CoS (live):  {len(gdf):,}")
print(f"  Ryde:        {len(ryde_near):,} of {len(ryde_gdf):,}")
for src in rbg_gdf["source"].unique():
    n_near = len(rbg_near[rbg_near["source"] == src])
    n_total = len(rbg_gdf[rbg_gdf["source"] == src])
    print(f"  {src}: {n_near:,} of {n_total:,}")
print(f"  USYD:        {len(usyd_near):,} of {len(usyd_gdf):,}")
print(f"\nWithin {RADIUS_WIDE/1000:.0f} km of CBD (wide):")
print(f"  Ryde: {len(ryde_wide):,} of {len(ryde_gdf):,}")

# Combined GDF for the primary species map (species column preserved)
combined_near = pd.concat([
    gdf[["geometry", "source", "species"]],
    ryde_near[["geometry", "source", "species"]],
    rbg_near[["geometry", "source", "species"]],
    usyd_near[["geometry", "source", "species"]],
], ignore_index=True)
combined_near = geopandas.GeoDataFrame(combined_near, crs="EPSG:4326")
combined_near["species"] = [
    x if isinstance(x, str) else "unknown" for x in combined_near["species"]
]

# ---------------------------------------------------------------------------
# Map bounds and base geography
# ---------------------------------------------------------------------------

PAD = 0.005  # ~500 m edge padding

# CoS-tight bounds (maps 1, 3, 5, 6)
cos_minx, cos_miny, cos_maxx, cos_maxy = gdf.total_bounds
xlim = (cos_minx - PAD, cos_maxx + PAD)
ylim = (cos_miny - PAD, cos_maxy + PAD)

# Combined-data bounds for the primary species map (map 4)
all_minx, all_miny, all_maxx, all_maxy = combined_near.total_bounds
species_xlim = (all_minx - PAD, all_maxx + PAD)
species_ylim = (all_miny - PAD, all_maxy + PAD)

aus_poas = geopandas.read_file("geopandas-blog/aus_poas.shp")
# Suburb outlines for the CoS-tight crop
bbox_geom = box(xlim[0], ylim[0], xlim[1], ylim[1])
syd = aus_poas[aus_poas.intersects(bbox_geom)]
# Suburb outlines for the combined species map
syd_species = aus_poas[aus_poas.intersects(
    box(species_xlim[0], species_ylim[0], species_xlim[1], species_ylim[1])
)]

# ---------------------------------------------------------------------------
# Map 1 — all trees
# ---------------------------------------------------------------------------

print("\nRendering map: all trees…")
fig, ax = plt.subplots()
syd.plot(ax=ax, color="ghostwhite", edgecolor="black")
gdf.plot(ax=ax, color="darkgreen", marker=".", markersize=1, alpha=0.15)
clean_ax(ax)
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
clean_ax(ax)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
ax.set_title("Figs in the City of Sydney — live open data")
plt.tight_layout()
plt.savefig("out/figs_live.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved out/figs_live.png")

# ---------------------------------------------------------------------------
# Map 4 — PRIMARY: all species, all sources, unique colour+marker per species
# ---------------------------------------------------------------------------

print("Rendering PRIMARY MAP: all species (CoS + RBG + Centennial + Ryde)…")
markers = [".", "o", "v", "^", "<", ">", "1", "2", "3", "4",
           "8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|", "_"]
tab = list(mcolors.TABLEAU_COLORS.keys())
base = list(mcolors.BASE_COLORS.keys())
tab.extend(base)
colours = [c for c in tab if c != "w"]
idents = [{"m": m, "c": c} for m in markers for c in colours]
seed(42)
shuffle(idents)

tree_counts = combined_near["species"].value_counts()

fig_all, ax = plt.subplots()
syd_species.plot(ax=ax, color="ghostwhite", edgecolor="#cccccc", linewidth=0.5)
clean_ax(ax)
ax.set_xlim(*species_xlim)
ax.set_ylim(*species_ylim)

legend_handles = []
for i, (tree, count) in enumerate(tree_counts.items()):
    ident = idents[i % len(idents)]
    sub = combined_near[combined_near["species"] == tree]
    try:
        common = "—" + tree_common_names[tree][0].title()
    except (KeyError, IndexError):
        common = ""
    label = rf"$\it{{{tree}}}${common} ({count})" if count >= 50 else None
    sub.plot(ax=ax, color=ident["c"], marker=ident["m"],
             alpha=0.4, markersize=3)
    if label:
        legend_handles.append(
            mlines.Line2D([], [], color=ident["c"], marker=ident["m"],
                          linestyle="None", markersize=5,
                          label=rf"$\it{{{tree}}}${common} ({count})")
        )

ax.legend(handles=legend_handles, loc="upper left",
          bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
          prop={"size": 5}, labelspacing=0.5)
title = (
    f"{len(combined_near):,} Trees — City of Sydney, Royal Botanic Garden, "
    f"Centennial Parklands, City of Ryde (partial), University of Sydney"
)
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
    clean_ax(ax)
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
# Map 6 — tight CoS view showing where supplementary data sits relative to CoS
# ---------------------------------------------------------------------------

print("\nRendering map: CoS + nearby RBG/Centennial (CoS-tight crop)…")

# Combine to a slim GDF with just geometry + source
combined_tight = pd.concat([
    gdf[["geometry", "source"]],
    ryde_near[["geometry", "source"]],
    rbg_near[["geometry", "source"]],
    usyd_near[["geometry", "source"]],
], ignore_index=True)
combined_tight = geopandas.GeoDataFrame(combined_tight, crs="EPSG:4326")

source_styles = {
    "City of Sydney (live)": {"color": "darkgreen", "ms": 1,  "alpha": 0.2,  "z": 2},
    "Centennial Parklands":  {"color": "C2",        "ms": 2,  "alpha": 0.55, "z": 3},
    "Royal Botanic Garden":  {"color": "C0",        "ms": 3,  "alpha": 0.7,  "z": 4},
    "City of Ryde (2013)":   {"color": "C1",        "ms": 1,  "alpha": 0.3,  "z": 2},
    "Mt Annan":              {"color": "C4",        "ms": 3,  "alpha": 0.7,  "z": 3},
    "Mt Tomah":              {"color": "C5",        "ms": 3,  "alpha": 0.7,  "z": 3},
    "University of Sydney":  {"color": "C3",        "ms": 2,  "alpha": 0.6,  "z": 3},
}

fig_t, ax = plt.subplots()
syd.plot(ax=ax, color="ghostwhite", edgecolor="black", linewidth=0.5)
legend_handles = []
for src, style in source_styles.items():
    sub = combined_tight[combined_tight["source"] == src]
    if not len(sub):
        continue
    sub.plot(ax=ax, color=style["color"], marker=".",
             markersize=style["ms"], alpha=style["alpha"], zorder=style["z"])
    legend_handles.append(mpatches.Patch(color=style["color"], label=f"{src} ({len(sub):,})"))
ax.legend(handles=legend_handles, loc="lower right", fontsize=7)
clean_ax(ax)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
ax.set_title(f"Trees near Sydney CBD — CoS view ({len(combined_tight):,} within crop)")
plt.tight_layout()
plt.savefig("out/combined_cos_crop.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved out/combined_cos_crop.png")

# ---------------------------------------------------------------------------
# Map 7 — 10 km wide view, all sources
# ---------------------------------------------------------------------------

print("Rendering map: combined data — 10 km view…")

combined_wide = pd.concat([
    gdf[["geometry", "source"]],
    ryde_wide[["geometry", "source"]],
    rbg_near[["geometry", "source"]],
    usyd_near[["geometry", "source"]],
], ignore_index=True)
combined_wide = geopandas.GeoDataFrame(combined_wide, crs="EPSG:4326")

# Derive map extent from the 20 km projected buffer
cbd_buffer_4326 = (
    geopandas.GeoDataFrame(geometry=[cbd_proj.buffer(RADIUS_WIDE + 500)], crs="EPSG:28356")
    .to_crs(epsg=4326)
)
bx0, by0, bx1, by1 = cbd_buffer_4326.total_bounds
wide_xlim = (bx0, bx1)
wide_ylim = (by0, by1)

wide_bbox = box(wide_xlim[0], wide_ylim[0], wide_xlim[1], wide_ylim[1])
syd_wide = aus_poas[aus_poas.intersects(wide_bbox)]

fig_w, ax = plt.subplots()
syd_wide.plot(ax=ax, color="ghostwhite", edgecolor="black", linewidth=0.5)
legend_handles = []
for src, style in source_styles.items():
    sub = combined_wide[combined_wide["source"] == src]
    if not len(sub):
        continue
    sub.plot(ax=ax, color=style["color"], marker=".",
             markersize=style["ms"], alpha=style["alpha"], zorder=style["z"])
    legend_handles.append(mpatches.Patch(color=style["color"], label=f"{src} ({len(sub):,})"))
ax.legend(handles=legend_handles, loc="lower right", fontsize=7)
clean_ax(ax)
ax.set_xlim(*wide_xlim)
ax.set_ylim(*wide_ylim)
ax.set_title(f"Trees within 10 km of Sydney CBD — {len(combined_wide):,} total")
plt.tight_layout()
plt.savefig("out/combined_10km.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved out/combined_10km.png")

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

print("\nAll done.  Outputs in out/")
print("\nLive dataset column summary:")
print(gdf.describe(include="all").T[["count", "unique", "top"]].to_string())
