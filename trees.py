#%%
"""
Plot trees in sydney on a map
Data from here: https://www.righttoknow.org.au/request/tree_data#incoming-6569
"""
#%%
from datetime import datetime
from matplotlib import font_manager as fm, rcParams
from pandas.api.types import CategoricalDtype
from random import shuffle
from shapely.geometry import Point, LineString, Polygon
from tree_common_to_latin_names_map import tree_common_names
import geopandas
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
import requests
import shapely
import time

#%%
res_multiplier = 1
cm_to_inch = 2.54
plt.rcParams["figure.figsize"] = [
    (39 / cm_to_inch) * res_multiplier,
    (49 / cm_to_inch) * res_multiplier,
]

#%% This shows that we're realing with coords and a species.
# There isn't a tree ID or any other tree metadata, but that might be a feature
# of the wording of the FOI request that produced the data.
# There's a LOT of trees!.
df = pd.read_excel(os.path.join("in", "Street tree data.xlsx"))
df.feat_cent_east = df.feat_cent_east.astype(float)
df.feat_cent_north = df.feat_cent_north.astype(float)
print(df.columns)
print(df.shape)
df.sample(5)
#%%
print(f"There are {len(df.species.unique())} species of tree covered by this data")
top = 30
df.species.value_counts()[:top].plot(rot=90, kind="bar")
plt.title(f"Top {top} most common trees in Sydney")
plt.savefig(f"Top {top} most common trees in Sydney", bbox_inches="tight")
#%%
loners = df.species.value_counts()[df.species.value_counts() == 1]
print(len(loners), ", ".join(loners.index))

#%% These numbers are in northings and eastings,
# They need to be converted to lat long before any of the normal ways of dealing
# with these things become easy. As they're in Sydney, they are covered by MGA56:
#     GDA94 / MGA zone 56
#     EPSG::28356
#     ProjectedCRS
#     Valid
#     Australia - onshore and offshore between 150°E and 156°E.
# Nice lat/long uses WGS84, aka epsg 4326.
# If you want to look this stuff up, this site is useful, if a bit clunky:
# http://www.epsg-registry.org/
# If we plot that onto a simple scatter, we see some pretty crazy outliers.
def make_pt(row):
    p = Point(row.feat_cent_east, row.feat_cent_north)
    return p


def load_ryde_data():
    gdf = geopandas.read_file(os.path.join("in", "ryde", "public-trees-2013.shp"))
    geometry = gdf.geometry
    geometry.crs = {"init": "epsg:28356", "no_defs": True}
    geometry = geometry.to_crs(epsg=4326)
    gdf["geometry"] = geometry
    gdf["species"] = gdf.apply(lambda x: "unknown")
    # gdf.plot(column="Height", legend=True, markersize=2)
    return gdf


world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
aus = world[world.name == "Australia"].plot()

geometry = geopandas.GeoSeries(df.apply(make_pt, axis=1))
geometry.crs = {"init": "epsg:28356", "no_defs": True}
geometry = geometry.to_crs(epsg=4326)
gdf = geopandas.GeoDataFrame(df)
gdf["geometry"] = geometry
gdf = gdf.append(load_ryde_data())
gdf.plot(ax=aus, c="r")

#%% Let's work out where we are in the world.
# This data is supposed to cover Sydney, so if we make a radius around the
# centroid of the 2000 postcode we should be in the right ballpark.
aus_poas = geopandas.read_file("geopandas-blog/aus_poas.shp")
syd2000 = aus_poas.query("code == 2000")
middle = syd2000.centroid.iloc[0]
syd2000.plot()

#%% This histogram shows that almost all of the points are very close, with only
# a tiny number of rogues off in the middle of nowhere.
gdf["distance_from_middle"] = gdf.geometry.distance(middle)
print(gdf.distance_from_middle.max())
gdf.distance_from_middle.hist(bins=30)
#%% if we look at the outliers we'll see that 2 of them just have zeroes, and
# one has probably been typed in missing some numbers, probably should be
# 11333093.0, not 333093.0
gdf[gdf.distance_from_middle > 0.3]
#%% But, as it's only a few trees, and we don't have tree IDs to be able to go
# back and find out what it's supposed to be, let's just can them.
cropped_gdf = gdf[gdf.distance_from_middle < 1]

#%% If we plot that onto a filtered set of suburbs:
aus_poas["distance_from_middle"] = aus_poas.geometry.distance(middle)
syd = aus_poas[aus_poas.distance_from_middle < 0.18]
syd_map = syd.plot(color="ghostwhite", edgecolor="black")
cropped_gdf.plot(ax=syd_map, color="r", marker=".", alpha=0.1)
plt.savefig("all_trees", bbox_inches="tight")


#%%
gdf.species = [x if type(x) is str else "unknown" for x in gdf.species]
#%% that's a lot of trees!
# What about just figs?
figs = {
    "Ficus microcarpa var hillii": {"colour": "C0", "marker": "+", "alpha": 0.3},
    "Ficus benjamina": {"colour": "C1", "marker": ".", "alpha": 0.3},
    "Ficus rubiginosa": {"colour": "C2", "marker": "1", "alpha": 0.3},
    "Ficus macrophylla": {"colour": "C3", "marker": "*", "alpha": 0.9},
    "Ficus sp.": {"colour": "C4", "marker": "p", "alpha": 0.3},
    "Ficus carica": {"colour": "C5", "marker": "o", "alpha": 0.3},
    "Ficus elastica": {"colour": "C6", "marker": "s", "alpha": 0.3},
}

fig_df = gdf[["Ficus" in x for x in gdf.species]]
fig_df["colour"] = fig_df.apply(lambda row: figs[row.species]["colour"], axis=1)
fig_df["alpha"] = fig_df.apply(lambda row: figs[row.species]["alpha"], axis=1)
fig_df["marker"] = fig_df.apply(lambda row: figs[row.species]["marker"], axis=1)
print(f"There are {fig_df.shape[0]} figs in this area.")
fig_df.head()

#%%
syd_map = syd.plot(color="ghostwhite", edgecolor="black")
fig_df.plot(ax=syd_map, color=fig_df.colour, marker=".", alpha=0.6)
# list(gdf.marker), alpha=list(gdf.alpha)) # it hates me :(
plt.legend(
    handles=[mpatches.Patch(color=figs[f]["colour"], label=f) for f in figs.keys()],
    loc="upper left",
)
plt.title("Figs in Sydney (as of October 17, 2016)")
plt.savefig("fig_tree_map", bbox_inches="tight")


#%% can we make that red splodge of trees any easier to understand?
tree_counts = gdf.species.value_counts()
markers = [  # see https://matplotlib.org/api/markers_api.html
    ".",
    ",",
    "o",
    "v",
    "^",
    "<",
    ">",
    "1",
    "2",
    "3",
    "4",
    "8",
    "s",
    "p",
    "P",
    "*",
    "h",
    "H",
    "+",
    "x",
    "X",
    "D",
    "d",
    "|",
    "_",
]
tab = list(mcolors.TABLEAU_COLORS.keys())
base = list(mcolors.BASE_COLORS.keys())
tab.extend(base)
colours = tab
colours.remove("w")
idents = []
for m in markers:
    for c in colours:
        idents.append({"m": m, "c": c})
shuffle(idents)
print(f"there are {len(idents)} possible idents and {len(tree_counts)} types of tree")
#%%
syd_map = syd.plot(color="ghostwhite", edgecolor="white")
plot_tight_on_CoS = False
if plot_tight_on_CoS:
    plt.xlim((151.17, 151.28))
    plt.ylim((-33.925, -33.845))
else:
    plt.xlim((151.05, 151.3))
    plt.ylim((-33.95, -33.75))
# print(plt.axis())
for i, x in enumerate(tree_counts.iteritems()):
    tree = x[0]
    count = x[1]
    temp_df = gdf[gdf.species == tree]
    try:
        common_name = "—" + tree_common_names[tree][0].title()
    except:
        common_name = ""
    label = r"$\it{" + tree + "}$" + f"{common_name} ({count})" if count >= 50 else None
    temp_df.plot(
        ax=syd_map,
        color=idents[i]["c"],
        marker=idents[i]["m"],
        alpha=0.4,
        label=label,
        markersize=3,
    )

plt.legend(loc="upper right", markerscale=2, prop={"size": 6}, labelspacing=0)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
title = f"{gdf.shape[0]} Trees in the City of Sydney (2017 data)"
plt.title(title)
metadata = {
    "Title": title,
    "Author": "Ben Doherty",
    "Description": """All trees in the City of sydney, as derived from a freedom 
                            of information act request in 2017 by Luke Bacon (@equivalentideas)
                            https://www.righttoknow.org.au/request/tree_data#incoming-6569
                            
                            It misses the trees managed by Parks NSW, USYD, and the 
                            Royal Botanic Gardens, as well as trees on private property.""",
    "Copyright": "",
    "Creation Time": "Fri, 7 Feb 2020 23:45:00 AEST",
    "Software": "Matplotlib",
    "Disclaimer": "",
    "Warning": "",
    "Source": "Dell Inspiron, VS Code",
    "Comment": "",
}
plt.savefig(
    "all_tree_map.png", bbox_inches="tight", dpi=600, metadata=metadata, format="png"
)
plt.savefig(
    "all_tree_map.svg", bbox_inches="tight", dpi=600, metadata=metadata, format="svg"
)

#%%
# import contextily as ctx
# tdf = geopandas.read_file(geopandas.datasets.get_path('nybb'))
# ax = tdf.plot(figsize=(10, 10), alpha=0.5, edgecolor='k')
# ctx.add_basemap(ax, url=ctx.providers.Stamen.TonerLite)
# ax.set_axis_off()
# %%
# this data is from the opentrees.org website
SB_df = pd.read_csv(os.path.join("in", "alltrees.csv"))
SB_df.sample(10)  # .to_html()
