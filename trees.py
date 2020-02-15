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
from shapely.geometry import Point
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
df = pd.read_excel("Street tree data.xlsx")
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
loners = df.species.value_counts()[df.species.value_counts()==1]
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


world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
aus = world[world.name == "Australia"].plot()

geometry = geopandas.GeoSeries(df.apply(make_pt, axis=1))
geometry.crs = {"init": "epsg:28356", "no_defs": True}
geometry = geometry.to_crs(epsg=4326)
gdf = geopandas.GeoDataFrame(df)
gdf["geometry"] = geometry
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
gdf[gdf.distance_from_middle > 0.1]
#%% But, as it's only a few trees, and we don't have tree IDs to be able to go
# back and find out what it's supposed to be, let's just can them.
gdf = gdf[gdf.distance_from_middle < 0.1]

#%% If we plot that onto a filtered set of suburbs:
aus_poas["distance_from_middle"] = aus_poas.geometry.distance(middle)
syd = aus_poas[aus_poas.distance_from_middle < 0.05]
syd_map = syd.plot(color="ghostwhite", edgecolor="black")
gdf.plot(ax=syd_map, color="r", marker=".", alpha=0.1)
plt.savefig("all_trees", bbox_inches="tight")

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

#%%
tree_common_names = {
"Lophostemon confertus":
    ["vinegartree", "Brisbane boxtree", "Brisbane-box", "Queensland-box", "vinegar-tree", "brushbox", "red-box", "vinegartree", "Brisbane boxtree", "Brisbane-box", "Queensland-box", "vinegar-tree", "brushbox", "red-box"],

"Platanus acerifolia":
    ["London planetree", "London planetree"],

"Melaleuca quinquenervia":
    ["Japanese paper wasp", "Mao-Holzrose", "aceite de cayeput", "ahambo", "balsamo de cayeput", "belbowrie", "bottle brush tree", "broad-leaved paperbark tree", "broadleaf paperbark tree", "broadleaf teatree", "cajeput", "capeputi", "corcho", "five-veined paperbark tree", "itahou", "kayu putih", "kinindrano", "melaleuca", "niaouli", "niaouli", "numbah", "oli", "paper bark tree", "paperbark teatree", "punk tree", "white bottlebrush tree", "bottle brush tree", "cajeput tree", "melaleuca", "niaouli", "paperbark", "punktree"],

"Tristaniopsis laurina":
    ["water-gum", "kanuka", "water-gum", "water-gum", "water-gum", "kanuka"],

"Robinia pseudoacacia Frisia":
    ["Golden Robinia"],

"Jacaranda mimosifolia":
    ["Black poui", "Flamboyán azul", "Jacaranda"],

"Corymbia maculata":
    ["spotted gum", "spotted gum"],

"Cupaniopsis anacardioides":
    ["carrot weed", "carrotwood", "tuckeroo", "carrotwood", "cashew-leaf cupania", "beach-tamarind", "green-leaf-tamarind", "Carrotwood", "carrot weed", "carrotwood", "tuckeroo", "carrot weed", "carrotwood", "carrotwood", "cashew-leaf cupania", "beach-tamarind", "green-leaf-tamarind", "tuckeroo", "carrot weed", "carrotwood", "Carrotwood"],

"Elaeocarpus reticulatus":
    ["blue oliveberry", "fringetree", "lily-of-the-valley-tree", "scrub-ash", "ash quandong", "fairy petticoats", "blueberry-ash", "blue oliveberry", "fringetree", "lily-of-the-valley-tree", "scrub-ash", "ash quandong", "fairy petticoats", "blueberry-ash"],

"Callistemon viminalis":
    ["weeping bottlebrush", "drooping bottlebrush", "red bottlebrush", "creek bottlebrush", "weeping bottlebrush", "Limpiatubos lloron", "Rutenförmiger Zylinderputzer", "Limpiatubos lloron", "Rutenförmiger Zylinderputzer"],

"Ficus microcarpa var hillii":
    ["Hill's Weeping Fig"],

"Pistacia chinensis":
    ["Chinese pistache",  "Chinese pistache", "Chinese pistachio", "agiao", "huang lian mu", "Chinese pistache", "Chinese pistachio", "sanguido", "Chinese pistache", "Chinese pistachio", "Chinese pistache", "Chinese pistache"],

"Fraxinus griffithii":
    ["Formosan ash", "Griffith's ash", "Philippine ash", ";", "Formosan ash", "guang la shu", "Griffith's ash", "Philippine ash", ";"],

"Sapium sebiferum":
    ["Chinese tallow tree", "Florida aspen", "chicken tree", "popcorn tree", "tallowtree", "vegetable tallow", "white wax berry", "Chinese tallow tree", "Florida aspen", "chicken tree", "popcorn tree", "tallowtree", "vegetable tallow", "white wax berry"],

"Liquidambar styraciflua":
    ["Sweetgum", "sweetgum", "Sweetgum", "sweetgum", "American sweet gum", "sweet gum", "alligator-wood", "satin-walnut", "American-storax", "red-gum", "sweet-gum"],

"Lagerstroemia indica":
    ["Crepe-myrtle"],

"Eucalyptus microcorys":
    ["Australian tallowwood", "tallowwood", "Australian tallowwood", "tallowwood", "tallowwood", "eucalipto"],

"Flindersia australis":
    ["Australian-teak", "crow's-ash", "Australian-teak", "crow's-ash"],

"Celtis australis":
    ["European hackberry", "European nettletree", "Mediterranean hackberry", "honey-berry"],

"Populus simonii":
    ["Simon poplar",  "Simon poplar", "Chinese cottonwood", "Simon poplar", "Simon poplar"],

"Koelreuteria paniculata":
    ["goldenrain-tree", "varnishtree", "pride-of-India", "Blasenesche", "Goldenrain tree"],

"Fraxinus pennsylvanica":
    ["green ash",  "Green ash, Red ash", "green ash", "green ash", "red ash", "downy ash", "green ash", "downy ash", "green ash", "northern red ash", "red ash"],

"Angophora costata":
    ["smooth-bark-apple", "smooth-bark-apple"],

"Platanus orientalis Digitata":
    ["Cut-leaf Plane", "Cyprian Plane"],

"Eucalyptus sideroxylon":
    ["red ironbark", "red ironbark eucalyptus", "mugga ironbark", "red ironbark", "mugga ironbark", "three-fruit red ironbark", "black ironbark", "Red ironbark", "Red ironbark", "red ironbark"],

"Ulmus parvifolia":
    ["Chinese elm", "lacebark", "lacebark elm", "Chinese elm, Lacebark elm",   "Chinese elm", ";", "Chinese elm", "lacebark elm", "Chinese elm", "lacebark", "leatherleaf elm", "Chinese elm", "chamneureupnamu", "lang yu", "aki-nire", "lacebark", "lacebark elm", "Chinese elm", "lacebark elm", ";", "Chinese elm", "Chinese elm", "lacebark elm"],

"Livistona australis":
    ["Australian cabbage palm", "Australian fan palm", "Gippsland palm", "cabbage fan palm", "cabbage palm", "cabbage-tree palm", "cabbagetree", "Australian cabbage palm", "Australian fan palm", "Gippsland palm", "cabbage fan palm", "cabbage palm", "cabbage-tree palm", "cabbagetree"],

"Corymbia citriodora":
    ["lemonscented gum"],

"Magnolia grandiflora":
    [ "bull-bay", "laurier tulipier", "southern magnolia", "southern magnolia", "taisan-boku", "bull bay", "southern magnolia", "bull bay", "southern magnolia"],

"Populus nigra Italica":
    ["Lombardy poplar", "Lombardy Poplar", "Pyramidenpappel"],

"Eucalyptus sp.":
    [],

"Eucalyptus botryoides":
    ["southern mahogany", "Gippsland-mahogany", "southern-mahogany", "bangalay", "southern mahogany", "Gippsland-mahogany", "southern-mahogany", "bangalay", "pu tao an"],

"Fraxinus oxycarpa Raywoodii":
    ["Claret Ash"],

"Ficus benjamina":
    ["weeping fig", "benjamin fig", "ficus tree"],

"Eucalyptus scoparia":
    ["wallangarra white gum", "Wallangarra white gum", "willow gum", "Wallangarra white gum", "willow gum"],

"Liriodendron tulipifera":
    ["tulip tree", "American tulip tree", "tulipwood", "tuliptree", "tulip poplar", "whitewood", "fiddletree", "yellow-poplar"],

"Ficus rubiginosa":
    ["Port Jackson fig", "little leaf fig", "rusty fig", "rusty-leaved fig", "Higuera mohosa", "Port Jackson fig", "Rostrote Feige", "Higuera mohosa", "Port Jackson fig", "Rostrote Feige", "Port Jackson fig", "little leaf fig", "rusty fig", "rusty-leaved fig", "Port Jackson fig", "Port Jackson fig", "littleleaf fig", "Illawarra fig", "rusty fig", "Port Jackson fig", "larger small-leaf fig", "small-leaf fig", "littleleaf fig", "Illawarra fig", "rusty fig", "Port Jackson fig", "larger small-leaf fig", "small-leaf fig", "littleleaf fig"],

"Olea europaea subsp europaea":
    ["European olive", "European olive", "olive", "mu xi lian", "olive"],

"Celtis occidentalis":
    ["common hackberry", "hackberry", "western hackberry"],

"Platanus orientalis":
    [ "Oriental Plane", "Orientalische Platane", "Oriental planetree", "Oriental plane", "Oriental planetree", "chenar", "san qiu xuan ling mu", "Oriental plane", "Oriental planetree", "chenar", "Oriental planetree"],

"Waterhousea floribunda":
    ["weeping lilly pilly", "weeping myrtle", "weeping lilly pilly", "weeping myrtle"],

"Callistemon salignus":
    ["stonewood", "white bottlebrush", "willow bottlebrush", "white bottlebrush", "willow bottlebrush", "stonewood", "pink-tip bottlebrush", "pink-tips"],

"Brachychiton acerifolius":
    ["flame bottletree", "Illawara flametree", "flame bottletree", "flame kurrajong", "lace-bark tree", "Australian flametree", "Illawarra flametree", "flame bottletree", "flame kurrajong", "flametree", "white crowsfoot", "lacebarktree", "flame bottletree", "Australian flametree", "Illawarra flametree", "flame bottletree", "flame kurrajong", "flametree", "white crowsfoot", "lacebarktree"],

"Syzygium paniculatum":
    ["magenta lilly-pilly", "brush cherry", "Australian brush-cherry", "magenta lilly-pilly", "brush cherry"],

"Harpephyllum caffrum":
    ["Kaffir Plum", "Wild Plum", "Kaffir Plum", "Wild Plum", "Kaffir-date", "Kaffir-plum", "Kaffir-date", "Kaffir-plum"],

"Casuarina cunninghamiana":
    ["Cunningham's beefwood", "river she-oak", "beefwood", "pinheiro-australiano", "river she-oak", "river-oak", "casuarina", "creek-oak", "fire-oak", "Cunningham's beefwood", "river she-oak"],

"Melaleuca styphelioides":
    ["prickly-leaf teatree", "prickly-leaf teatree", "prickly paperbark", "prickly-leaf paperbark", "prickly-leaf teatree", "prickly paperbark", "prickly-leaf paperbark", "prickly-leaf teatree"],

"Populus deltoides":
    ["Eastern cottonwood"],

"Banksia integrifolia":
    ["coast banksia", "white-honeysuckle", "coast banksia", "white-honeysuckle", "white banksia", "white bottlebrush", "white-honeysuckle", "honeysuckle-oak", "coast banksia", "white banksia", "white bottlebrush", "white-honeysuckle", "honeysuckle-oak", "coast banksia"],

"Casuarina glauca":
    ["Brazilian beefwood", "gray she-oak", "longleaf ironwood", "saltmarsh ironwood", "scaly-bark beefwood", "suckering Australian-pine", "swamp oak", "swamp she-oak", "swamp oak", "Brazilian-oak", "longleaf ironwood", "marsh she-oak", "scaly-bark beefwood", "swamp she-oak", "swamp-oak", "saltmarsh ironwood", "grey bull-oak", "grey she-oak", "Gray sheoak"],

"Hymenosporum flavum":
    ["Queensland frangipani", "Queensland frangipani"],

"Corymbia eximia":
    ["yellow bloodwood", "yellow bloodwood", "yellow bloodwood", "yellow bloodwood", "yellow bloodwood"],

"Stenocarpus sinuatus":
    ["Queensland fire-wheel-tree", "tulip-flower", "tuliptree", "wheel-of-fire", "wheel-of-fire-tree", "wheeltree", "white beefwood", "white silky-oak", "fire-wheel-tree", "firetree", "white-oak", "firewheeltree", "firewheeltree", "Queensland fire-wheel-tree", "tulip-flower", "tuliptree", "wheel-of-fire", "wheel-of-fire-tree", "wheeltree", "white beefwood", "white silky-oak", "yiel-yiel", "fire-wheel-tree", "firetree", "white-oak"],

"Glochidion ferdinandii":
    ["cheese tree"],

"Acacia binervia":
    ["coast myall", "coast myall"],

"Agonis flexuosa":
    ["willow myrtle", "willow-peppermint", "Western Australian myrtle", "Western Australian peppermint", "Western Australian willow myrtle", "willow myrtle", "willow-peppermint", "Western Australian myrtle", "Western Australian peppermint", "Western Australian willow myrtle"],

"Waterhousea flori. 'Green Ave'":
    [],

"Syzygium luehmannii":
    ["riberry", "small leaved lilly pilly", "cherry satinash", "cherry alder",  "clove lilli pilli"],

"Backhousia citriodora":
    ["Australian lemon myrtle", "lemon ironwood", "lemon-scent backhousia", "lemon-scent myrtle", "lemon-scent verbena", "sweet verbena myrtle", "sweet verbena-tree", "Australian lemon myrtle", "lemon ironwood", "lemon-scent backhousia", "lemon-scent myrtle", "lemon-scent verbena", "sweet verbena myrtle", "sweet verbena-tree"],

"Gleditsia trican Sunburst":
    [],

"Eucalyptus saligna":
    ["Sydney bluegum", "blue gum", "Sydney blue gum", "liu ye an", "blue gum", "Sydney blue gum", "Sydney bluegum"],

"Castanospermum australe":
    ["Australian Chestnut", "Black Bean", "Moreton Bay Chestnut", "Australian chestnut", "Moreton Bay chestnut", "Australian Chestnut", "Black Bean", "Moreton Bay Chestnut", "Australian-chestnut", "Moreton Bay-bean", "Moreton Bay-chestnut", "beantree", "black-bean", "Australian-chestnut", "Moreton Bay-bean", "Moreton Bay-chestnut", "beantree", "black-bean"],

"Pittosporum rhombifolium":
    ["Queensland pittosporum", "white-holly", "diamond pittosporum", "diamond-leaf pittosporum", "diamond-leaf-laurel", "small-fruit pittosporum", "hollywood", "Queensland pittosporum", "white-holly", "diamond pittosporum", "diamond-leaf pittosporum", "diamond-leaf-laurel", "small-fruit pittosporum", "hollywood"],

"Pyrus calleryana":
    ["Callery pear", "dou li", "mame-nashi",  "callery pear", "Bradford pear", "Callery pear", "Callery pear",  "Callery pear", "dou li", "mame-nashi", "Callery pear", "Bradford pear", "callery pear", "dou li", "callery pear", "mame-nashi", "dou li", "Bradford pear", "Callery pear"],

"Schinus molle":
    ["California peppertree", "Peruvian peppertree", "Peruvian-mastictree", "peppertree"],

"Acer negundo":
    ["Ashleaf Maple", "ash-leaved maple, box-elder",  "Box-elder, Ash-leaved maple", "ashleaf maple", "California boxelder", "Manitoba maple", "ash-leaf maple", "box elder", "boxelder", "boxelder maple", "three-leaf maple", "western boxelder", "box elder", "fu ye feng", "klen âsenelistnyj", "ash-leaf maple", "box-elder", "tonerikoba-no-kaede", "three-leaf maple", "negundodanpung", "California box-elder", "Manitoba maple", "ash-leaved maple", "box-elder", "box-elder maple", "three-leaved maple", "western box-elder", "ash-leaved maple, box-elder"],

"Robinia pseudoacacia":
    ["Black locust", "False acacia", "False-acacia", "Post locust", "acacia blanc", "black locust", "false acacia", "robinia akacjowa", "robinier", "robinier faux acacia", "robinier faux-acacia", "yellow locust"],

"Acer Buergeranum":
    [],

"Acmena smithii":
    [],

"Melaleuca armillaris":
    ["bracelet honey myrtle", "giant honey myrtle", "Dropping Melaleuca", "Prickly-leaved tea-tree", "Dropping Melaleuca", "Prickly-leaved tea-tree"],

"Washingtonia robusta":
    ["Mexican fan palm", "Washington fan palm"],

"Callistemon citrinus":
    ["crimson bottlebrush", "Crimson bottlebrush", "Crimson bottlebrush", "calistemone", "crimson bottlebrush", "escova-de-garrafa"],

"Podocarpus elatus":
    ["Brown pine", "Illawarra plum", "Plum pine", "Yellow pine", "plum-pine", "pencil-cedar", "she-pine", "white-pine", "white-plum", "yellow-pine", "brown-pine", "plum-pine", "pencil-cedar", "she-pine", "white-pine", "white-plum", "yellow-pine", "brown-pine", "Brown pine", "Illawarra plum", "Plum pine", "Yellow pine"],

"Cinnamomum camphora":
    ["Japanese camphor", "alcanfor", "alcanforero", "arvore da camphora", "campher", "camphor laurel", "camphor tree", "camphre", "camphrier", "canfora", "kampferbaum", "kuso-no-ki", "camphor laurel", "camphor tree", "camphortree"],

"Casuarina sp.":
    [],

"Gleditsia trican Shademaster":
    [],

"Eucalyptus robusta":
    ["swamp messmate", "swamp stringybark", "swamp-mahogany", "Swampmahogany", "swampmahogany"],

"Eucalyptus nicholii":
    ["narrow-leaved peppermint", "willow peppermint", "narrow-leaf black peppermint", "narrow-leaf peppermint", "willow peppermint", "willow-leaf peppermint", "small-leaf peppermint", "narrow-leaf black peppermint", "narrow-leaf peppermint", "willow peppermint", "willow-leaf peppermint", "small-leaf peppermint"],

"Quercus palustris":
    ["pin oak", "pin oak", "pin oak", "swamp oak", "swamp pin oak", "pin oak"],

"Harpulia pendula":
    ["Moreton Bay tulipwood", "Queensland tulipwood", "black tulipwood", "tulip lancewood", "tulipwood", "black-tulip", "mogun-mogun", "Moreton Bay tulipwood", "Queensland tulipwood", "black tulipwood", "tulip lancewood", "tulipwood", "black-tulip", "mogun-mogun"],

"Brachychiton discolor":
    ["Queensland lacebark", "scrub bottletree", "lace kurrajong", "lacebarktree", "lacewoodtree", "brush kurrajong", "pink kurrajong", "white kurrajong", "white-poplar", "hat-tree", "Queensland lacebark", "scrub bottletree", "lace kurrajong", "lacebarktree", "lacewoodtree", "brush kurrajong", "pink kurrajong", "white kurrajong", "white-poplar", "hat-tree"],

"Fraxinus sp.":
    [],

"Cyathea cooperi":
    ["Australian tree fern", "Cooper's cyathea", "fanjan Australien", "fougère arborescente d'Australia", "lacy tree fern", "scaly tree fern", "straw tree fern", "Australian tree fern", "Cooper's cyathea", "fanjan Australien", "fougère arborescente d'Australia", "lacy tree fern", "scaly tree fern", "straw tree fern", "Cooper's cyathea", "scaly tree fern", "Australian tree fern", "highland lace", "straw tree fern", "Australian tree fern", "highland lace", "straw tree fern", "Cooper's cyathea", "scaly tree fern"],

"Koelreuteria elegans subsp For":
    [],

"Grevillea robusta":
    ["chêne d'Australie", "grevillaire", "silk oak", "silky oak", "silver oak", "Roble australiano", "Roble de seda", "Silk oak", "chêne d'Australie", "grevillaire", "silk oak", "silky oak", "silver oak", "silkoak", "silky oak", "silver oak", "Australian silky-oak", "silk-oak", "silky-oak", "southern silky-oak", "shinobu-no-ki"],

"Melaleuca linariifolia":
    ["Cajeput tree", "Cajeput tree", "cajeput tree", "narrow-leaf paperbark", "snow-in-summer", "flax-leaf paperbark", "narrow-leaf teatree", "cajeput tree", "narrow-leaf paperbark", "snow-in-summer", "flax-leaf paperbark", "narrow-leaf teatree"],

"Plumeria acutifolia":
    [],

"Olea africana":
    [],

"Melia azedarach":
    ["Indian lilac", "Persian lilac", "Sichuan pagoda-tree", "alelaila", "amargoseira-do-Himalaio", "arbre à chapelets", "bakain", "chinaberry", "chuan liang zi", "dake", "indischer Zedrachbaum", "jazmin", "lelah", "lilas", "lilas de Perse", "lilas de l'Inde", "lilas des Indes", "margosa tree", "margosier", "melia", "paraíso", "para‘isu", "persischer Flieder", "petit lilas", "prais", "pride-of-India", "sendan", "sili", "sita", "syringa berrytree", "tili", "tira", "umbrella tree", "white cedar", "‘ilinia", "‘inia"],

"Murraya paniculata":
    ["Chinese box", "orange-jessamine", "Café de la India", "China Box", "Kamini", "Orange jessamine"],

"Leptospermum sp.":
    [],

"Ginkgo biloba":
    [ "ginkgo", "maidenhair-tree", "common ginkgo", "maidenhair tree"],

"Phoenix canariensis":
    ["Canary Island date palm", "Canary date palm", "dattier des Canaries", "phoenix palm",  "tamareira-das-Canárias", "Canary Island date palm", "Canary Island palm", "Canary date palm"],

"Syzygium australe":
    ["brush-cherry", "Australian brush-cherry", "Australian water-pear", "scrub-cherry", "creek lilly-pilly", "creek satin-ash", "creek-cherry", "brush-cherry", "Australian brush-cherry", "Australian water-pear", "scrub-cherry", "creek lilly-pilly", "creek satin-ash", "creek-cherry", "brush-cherry"],

"Syagrus romanzoffianum":
    ["queen palm", "giriba palm", "queen palm", "queen palm", "giriba palm", "queen palm"],

"Platanus occidentalis":
    ["American sycamore", "buttonball", "buttontree", "buttonwood", "eastern sycamore", "American plane", "American planetree", "American sycamore", "sycamore", "American sycamore", "sycamore",  "American sycamore", "American sycamore", "sycamore", "American sycamore", "eastern plane-tree", "eastern sycamore", "xuan ling mu", "buttonball", "buttontree", "buttonwood", "eastern sycamore", "American plane", "American planetree", "American sycamore", "sycamore", "sycamore", "American plane-tree", "American sycamore", "buttonball", "buttonball tree", "buttonwood"],

"Callistemon sp.":
    [],

"Metrosideros excelsa":
    ["New Zealand Christmas tree", "New Zealand Christmastree", "pohutukawa", "New Zealand Christmastree"],

"Archontophoenix cunninghamiana":
    ["piccabeen bangalow palm", "piccabeen palm", "bangalow palm", "Picabeen palm", "piccabeen bangalow palm", "piccabeen palm", "bangalow palm", "Picabeen palm"],

"Archontophoenix alexandrae":
    ["Alexandra palm", "northern bangalow palm", "king palm", "Picabeen palm", "Alexandra palm", "northern bangalow palm", "king palm", "Picabeen palm", "Alexandra palm", "Alexandra palm"],

"Lagunaria patersonii":
    ["cow itch tree", "cow itch tree"],

"Populus x canadensis":
    [],

"Hibiscus tiliaceus":
    [""],

"Prunus sp.":
    [],

"Eucalyptus grandis":
    ["grand eucalyptus", "rose gum", "saligna gum", "scrub gum", "rose gum", "da an", "flooded gum", "saligna gum", "scrub gum", "rose gum", "flooded gum"],

"Agathis robusta":
    ["Queensland kauri pine", "Smooth-barked kauri", "asong", "muwaka", "ogapa", "Queensland kauri", "Dundathu-pine", "Queensland kauri", "Queensland kauri-pine", "South Queensland kauri", "South Queensland kauri-pine", "kauri-pine", "smoothbark kauri", "North Queensland kauri", "North Queensland kauri-pine", "Dundathu-pine", "Queensland kauri", "Queensland kauri-pine", "South Queensland kauri", "South Queensland kauri-pine", "kauri-pine", "smoothbark kauri", "North Queensland kauri", "North Queensland kauri-pine", "Queensland kauri"],

"Fraxinus angustifolia":
    ["desert ash", "narrow-leaf ash", "desert ash", "narrow-leaf ash"],

"Grevillea sp.":
    [],

"Tristaniopsis 'Luscious'":
    [],

"Unknown":
    [],

"Araucaria cunninghamii":
    ["Hoop pine", "Queensland pine", "coonam", "coorong", "cumbertu", "Moreton Bay pine", "dorrigo pine", "Moreton Bay-pine", "Queensland-pine", "pinheiro-colonial", "pinheiro-cunnigami", "pinheiro-de-arco", "rocket-tree", "colonial-pine", "Richmond River-pine", "Dorrigo-pine", "hoop-pine", "nan yang shan"],

"Syzygium sp.":
    [],

"Corymbia ficifolia":
    ["redflower gum", "redflower gum", "red-flowering gum", "red-flower gum", "red-flower gum"],

"Macadamia integrifolia":
    ["Macadamia Nut", "Queensland Nut", "macadamia nut", "macadamia nut", "Queenslandnut", "ao zhou jian guo", "smooth-shell Queenslandnut", "macadamia-nut", "bopplenut", "poppelnut", "nut-oak", "Bauplenut", "bushnut", "macadamia nut", "Macadamia Nut", "Queensland Nut", "Queenslandnut", "smooth-shell Queenslandnut", "macadamia-nut", "bopplenut", "poppelnut", "nut-oak", "Bauplenut", "bushnut"],

"Persea americana":
    [],

"Afrocarpus Falcatus":
    ["bastard yellowwood", "outeniqua yellowwood", "common yellowwood", "Bastard yellowwood", "Outeniqua yellowwood", "inkoba", "mse mawe", "umGeya", "yellowwood", "African fern pine", "bastard yellowwood", "outeniqua yellowwood", "bastard yellowwood", "outeniqua yellowwood", "common yellowwood", "yellowwood", "Bastard yellowwood", "Outeniqua yellowwood", "inkoba", "mse mawe", "umGeya"],

"Zelkova serrata 'Green Vase'":
    [],

"Schefflera actinophylla":
    ["arbre ombrelle", "arbre pieuvre", "brassaia", "ivy palm", "octopus tree", "schefflera", "umbrella tree", "octopus tree", "umbrella tree"],

"Pyrus ussuriensis":
    [", ", "Ussurian pear", ", ", "Chinese pear", "Harbin pear", "Manchurian pear", "Ussuri pear", "Ussurian pear", "sandolbae", "qiu zi li", "Chinese pear", "Harbin pear", "Manchurian pear", "Ussuri pear", "Ussurian pear"],

"Pittosporum sp.":
    [],

"Leptospermum citratum":
    [],

"Fraxinus excelsior":
    ["European ash", "ash", "black ash", "ash"],

"Brachychiton sp.":
    [],

"Eucalyptus crebra":
    ["narrow-leaf ironbark", "red ironbark", "grey ironbark", "white ironbark", "narrow-leaf red ironbark", "narrowleaf red ironbark", "narrow-leaf ironbark", "red ironbark", "grey ironbark", "white ironbark", "narrow-leaf red ironbark", "narrowleaf red ironbark"],

"Quercus ilex":
    ["Evergreen Oak", "evergreen oak", "holly oak", "holm oak", "holly oak"],

"Celtis sinensis":
    ["Weeping Chinese hackberry",    "Chinese hackberry", "Japanese hackberry", "Chinese hackberry", "Japanese hackberry", "po shu", "enoki", "Chinese nettletree",  "Chinese hackberry", "Japanese hackberry", "Chinese nettletree"],

"Elaeocarpus eumundii":
    [],

"Acer sp.":
    [],

"Gordonia axillaris":
    ["Fried Egg Plant"],

"Prunus x blireana":
    ["purple-leafed plum" ,"double-flowering plum"],

"Citrus sp":
    [],

"Schefflera sp":
    [],

"Citharexylum spinosum":
    ["Florida fiddlewood", "fiddlewood", "masese", "spiny fiddlewood", "Florida fiddlewood", "fiddlewood", "masese", "spiny fiddlewood", "Florida fiddlewood", "fiddlewood", "spiny fiddlewood", "spiny fiddlewood", "fiddlewood", "Florida fiddlewood", "fiddlewood", "spiny fiddlewood", "spiny fiddlewood"],

"Cupressus sp.":
    [],

"Alectryon excelsum":
    ["titoki", "titoki"],

"Bizmarckia nobilis":
    [],

"Eucalyptus sideroxylon Rosea":
    [],

"Eucalyptus punctata":
    ["gray gum", "grey irongum", "long-cap grey gum", "grey gum", "grey iron gum", "long-cap grey gum", "grey gum", "grey iron gum", "gray gum"],

"Alectryon subcinereus":
    ["smooth rambutan", "bird's-eye", "wild quince", "smooth rambutan", "bird's-eye", "wild quince"],

"Bauhinia variegata":
    ["mountain-ebony", "orchidtree", "Bauhinia", "Camel's Foot", "Mountain Ebony", "Orchid tree", "Palo De Orquideas", "Pata De Rez", "Poor Man's Orchid", "St. Thomas' Tree", "Variegated Bauhinia",  "mountain ebony", "butterfly tree", "orchid tree", "purple orchid tree", "mountain-ebony", "orchidtree", "kachnar", "mountain ebony"],

"Camellia sp":
    [],

"Ficus macrophylla":
    ["Moreton Bay fig", "Australian banyan", "Moreton Bay fig", "permite", "black fig", "Australian banyan", "Moreton Bay fig", "black fig"],

"Populus simonii Fastigiata":
    [],

"Citrus limon":
    ["Lemon"],

"Hibiscus rosa sinensis":
    ["Chinese hibiscus", "China rose", "Hawaiian hibiscus", "rose mallow" , "shoeblackplant", ""],

"Melaleuca sp.":
    [],

"Brachychiton populneus":
    ["whiteflower kurrajong", "kurrajong", "kurrajong", "bottletree", "Whiteflower kurrajong", "Whiteflower kurrajong"],

"Podocarpus Falcatus":
    ["common yellowwood", "bastard yellowwood", "outeniqua yellowwood", "African fern pine", "weeping yew"],

"Leptospermum laevigatum":
    ["Australian myrtle", "Australian teatree", "coast teatree", "coastal teatree", "Australian teatree", "coastal teatree", "Australian teatree", "Australian myrtle", "coast teatree", "Australian myrtle", "Australian teatree", "coast teatree", "coastal teatree", "Australian teatree", "Australian teatree", "coastal teatree"],

"Pittosporum undulatum":
    ["Australian cheesewood", "Victorian box", "Victorian laurel", "mock orange", "native daphne", "orange pittosporum", "sweet pittosporum", "wild coffee", "Australian cheesewood", "Victorian box", "Victorian laurel", "mock orange", "sweet pittosporum"],

"Platanus insularis":
    ["Autumn Glory Plane"],

"Ulmus sp.":
    [],

"Acer platanoides":
    ["Norway maple", "Spitzahorn", "érable plane", "Norway maple", "Norway maple", "Norway Maple"],

"Alnus sp.":
    [],

"Populus alba":
    ["alamo blanco", "gattice", "gin-doro", "hakuyo", "peuplier blanc", "pioppo bianco", "silber-pappel", "silver-leaf poplar", "urajiro-hako-yanagi", "white poplar", "xin bai yang", "white poplar"],

"Hibiscus syriascus":
    [],

"Angophora floribunda":
    [],

"Leptospermum petersonii":
    ["common teatree", "common teatree", "lemon-scent teatree", "lemon-scent teatree"],

"Dypsis lutescens":
    ["yellow butterfly palm", "yellow butterfly palm", "areca palm", "areca-bambú", "cane palm", "golden-yellow palm", "palmeira-areca", "yellow palm", "yellow butterfly palm", "butterfly palm", "areca palm", "cane palm", "golden-yellow palm", "yellow palm", "yellow butterfly palm", "butterfly palm"],

"Camellia japonica":
    [ "Common camellia",   "Camellia", "Kamelie",  "common camellia", "camellia", "shan cha", "camellia", "camellia", "camellia", "Camellia", "Kamelie"],

"Eriobotrya japonica":
    ["Japanese medlar", "Japanese plum", "bibasse", "bibassier", "eriobotrya du Japon", "loquat", "néflier du Japon",  "Japanese medlar", "Japanische Mispel", "Loquat", "Nispero de Espana", "Wollmispel", "Japanese medlar", "Japanese plum", "bibasse", "bibassier", "eriobotrya du Japon", "loquat", "néflier du Japon", "loquat", "loquat", "loquat", "pi ba", "Japanese-medlar", "bibasse"],

"Howea forsteriana":
    ["sentrypalm", "Forster sentry palm", "palmeira-quência", "paradise palm", "sentry palm", "thatch-leaf palm", "kentia palm", "Forster sentry palm", "paradise palm", "sentry palm", "thatch-leaf palm", "kentia palm", "sentrypalm"],

"Banksia serrata":
    [],

"Nyssa sylvatica":
    ["black tupelo", "black-gum", "sour-gum", "black tupelo", "black-gum", "sour-gum", "black gum", "black tupelo", "blackgum", "black gum", "black tupelo", "pepperidge", "sour gum", "black gum", "black tupelo", "blackgum", "black tupelo", "black tupelo"],

"Photinia robusta":
    [],

"Mag. grandiflora 'Little Gem'":
    [],

"Chamaecyparis obtusa 'Crippsii":
    [],

"Celtis sp.":
    [],

"Platanus acerifolia 'Bloodgood":
    [],

"Acacia sp.":
    [],

"Olea europea var. africana":
    [],

"Eucalyptus globulus ssp. bicos":
    [],

"Eucalyptus gummifera":
    ["red bloodwood", "red bloodwood"],

"Murraya sp":
    [],

"Cupressus sempervirens":
    ["Italian cypress", "Mediterranean cypress", "mediteranski cempres", "obicni cempres"],

"Chamaecyparis lawsoniana":
    ["Lawson's-cypress", "ginger-pine", "Lawson's false cypress", "Port Orford-cedar", "Oregon-cedar", "Lawson cedar", "Oregon cedar", "Port Orford cedar", "Port Orford white cedar", "ginger-pine", "Oregon cedar", "Port Orford cedar", "Lawson's-cypress", "ginger-pine", "Lawson's false cypress", "mei guo bian bai", "Port Orford-cedar", "Oregon-cedar"],

"Albizia lophantha":
    [],

"Strelitzia nicolai":
    ["bird-of-paradise-tree", "white bird-of-paradise", "bird-of-paradise tree", "white bird-of-paradise tree", "bird-of-paradise tree", "white bird-of-paradise tree", "white bird-of-paradise", "bird-of-paradise-tree", "white bird-of-paradise"],

"Ficus sp.":
    [],

"Allocasuarina littoralis":
    ["black she-oak", "river black-oak", "bull-oak", "black she-oak", "river black-oak", "black she-oak", "river black-oak", "bull-oak", "black she-oak", "river black-oak"],

"Prunus cerasifera Nigra":
    [],

"Populus yunnanensis":
    ["Yunnan poplar", "dian yang", "Yunnan poplar"],

"Allocasuarina torulosa":
    ["forest-oak", "forest she-oak", "forest-oak", "river-oak", "rose she-oak", "forest-oak", "forest she-oak", "forest-oak", "river-oak", "rose she-oak"],

"Melaleuca bracteata":
    ["Oron-Wood", "Ridge Myrtle", "Oron-Wood", "Ridge Myrtle", "river teatree", "black teatree", "bracteate honey myrtle", "river teatree", "white cloudtree", "prickly-leaf teatree", "mock olive", "black teatree", "bracteate honey myrtle", "river teatree", "white cloudtree", "prickly-leaf teatree", "mock olive", "river teatree"],

"Washingtonia filifera":
    [ "California fan palm", "American cotton palm", "California washingtonia", "Washington palm", "cotton palm", "desert fan palm", "petticoat palm", "California fan palm", "California fan palm",  "California fan palm"],

"Photinia serrulata":
    [],

"Alectryon tomentosus":
    ["bed-jacket", "woolly rambutan", "red-jacket", "hairy alectryon", "hairy bird's-eye", "bed-jacket", "woolly rambutan", "red-jacket", "hairy alectryon", "hairy bird's-eye"],

"Howea belmoreana":
    ["Belmore sentrypalm", "Belmore sentrypalm", "Belmore palm", "Belmore sentry palm", "curly palm", "Belmore palm", "Belmore sentry palm", "curly palm"],

"Pittosporum eugenioides":
    ["lemonwood", "tarata lemonwood", "white mapau", "lemonwood", "tarata", "tarata lemonwood", "white mapau"],

"Mangifera sp.":
    [],

"Cotinous coggygria":
    [],

"Hakea sp.":
    [],

"Eucalyptus globulus":
    [],

"Alnus jorullensis":
    ["amieiro"],

"Catalpa bignonoides":
    [],

"Butia capitata":
    ["South American jelly palm", "pindo palm", "jelly palm", "South American jelly palm",   "South American jelly palm", "South American jelly palm", "pindo palm", "jelly palm"],

"Rose sp.":
    [],

"Araucaria heterophylla":
    ["Norfolk Island-pine", "Norfolk Island pine", "Star pine", "Norfolk Island pine", "Norfolk Island-pine", "Norfolk Island pine", "Star pine"],

"Allocasuarina cunninghamiana":
    [],

"Cordyline australis":
    ["cabbage tree", "cabbage tree", "New Zealand cabbagetree", "cabbage-palm", "cabbagetree", "fountain-dracaena", "giant-dracaena", "palm-lily", "grass-palm", "dracena-azul", "dracena-das-montanhas"],

"Palm sp.":
    [],

"Ulmus glabra Lutescens":
    [],

"Eucalyptus bicostata":
    [],

"Strelitzia sp.":
    [],

"Eucalyptus melliodora":
    ["yellow-box", "yellow ironbark", "yellow ironbox", "honey-box", "yellow-box", "yellow-box", "mi wei an", "yellow ironbark", "yellow ironbox", "honey-box", "yellow-box"],

"Eucalyptus ficifolia":
    ["redflower gum", "redflower gum"],

"Pyrus calleryana 'Chanticleer'":
    [],

"Phoenix roebelenii":
    ["pygmy date palm", "pygmy date palm", "Roebelin palm", "tamareira-anã", "tamareira-de-jardim", "miniature date palm", "pygmy date palm", "Roebelin palm", "miniature date palm", "pygmy date palm"],

"Araucaria columnaris":
    ["Cook's pine", "New Caledonia pine", "New Caledonia pine", "New Caledonia-pine", "Cook-pine", "New Caledonia-pine", "Cook-pine", "New Caledonia pine", "Cook's pine", "New Caledonia pine"],

"Acacia baileyana":
    ["Bailey's acacia", "Cootamundra wattle", "Bailey's acacia", "Cootamundra wattle", "Cootamundra wattle", "Bailey's acacia", "Bailey's wattle"],

"Eucalyptus mannifera":
    ["manna gum", "mountain spotted gum", "red spotted gum", "white brittle gum", "capertee brittle gum", "brittle gum", "broadleaf manna gum", "manna gum", "mountain spotted gum", "red spotted gum", "white brittle gum", "capertee brittle gum", "brittle gum", "broadleaf manna gum"],

"Nerium oleander":
    ["oleander", "rose-laurel", "rose bay"],

"Viburnum tinus":
    ["Laurustinus", "Immergrüner Schneeball", "Lorbeer-Schneeball", "Stein-Lorbeer", "laurustinus", "laurustinus", "Immergrüner Schneeball", "Lorbeer-Schneeball", "Stein-Lorbeer", "laurustinus", "laurustinus"],

"Eucalyptus saligna x botryoid":
    [],

"Populus alba Pyramidalis":
    [],

"Duranta erecta":
    ["golden dewdrops"],

"Ficus carica":
    ["Fig"],

"Liquidambar formosana":
    [ "Formosan gum", "Formosan-gum", "feng xiang shu", "Formosan sweet-gum",  "Formosan sweet-gum", "Formosan-gum"],

"Magnolia sp.":
    [],

"Populus sp.":
    [],

"Lophostemon suaveolens":
    ["paperbark-mahogany", "swamp terpentine", "swamp turpentine", "swamp-box", "swamp-mahogany", "water-gum", "paperbark-mahogany", "swamp terpentine", "swamp turpentine", "swamp-box", "swamp-mahogany", "water-gum"],

"Schefflera arboricola":
    ["Hawaiian-elf", "dwarf umbrella-tree", "miniature schefflera", "parasol-plant"],

"Trachycarpus fortunei":
    ["Chinese fan palm", "Chinese windmill palm", "chusan palm", "hemp palm", "hochstämmige Hanfpalme", "palma de jardín", "palmeira-moinho-de-vento-da-China", "palmier de Chine", "windmill palm", "Chinese windmill palm", "Chinese windmill palm", "windmill palm", "hemp palm", "Chusan fan palm", "Chusan palm", "Chinese fan palm", "Chinese windmill palm", "chusan palm", "hemp palm", "hochstämmige Hanfpalme", "palma de jardín", "palmeira-moinho-de-vento-da-China", "palmier de Chine", "windmill palm", "Chinese windmill palm", "windmill palm", "Chinese windmill palm", "palmeira-moinho-de-vento-da-China", "shuro", "windmill palm", "wa-juro", "hemp palm", "Chusan fan palm", "Chusan palm"],

"Ligustrum lucidum":
    ["broadleaf privet", "glossy privet", "large leaf privet", "ligustrum privet", "privet", "tree privet", "glossy privet"],

"Acer negundo Variegatum":
    [],

"Acer palmatum":
    ["Smooth Japanese maple",   "Japanese maple", "Japanese maple", "Japanese maple", "danpungnamu", "ji zhua feng", "takao-momiji", "smooth Japanese maple", "iroha-kaede", "iroha-momiji", "Japanese maple", "smooth Japanese maple", "Japanese maple"],

"Populus x canadensis Aurea":
    [],

"Syncarpia glomulifera":
    ["red luster", "turpentine", "turpentine-tree", "red turpentine", "turpentine tree", "turpentine tree", "red luster", "turpentine", "turpentine-tree", "red turpentine"],

"Eucalyptus tereticornis":
    ["Queensland blue gum", "blue gum", "forest red gum", "red gum", "red ironbark", "red irongum", "bastard-box", "mountain gum", "grey gum", "flooded gum", "stinking gum", "slaty gum", "forest redgum", "red ironbark", "Queensland blue gum", "blue gum", "forest red gum", "red gum", "red ironbark", "red irongum", "xi ye an", "bastard-box", "mountain gum", "grey gum", "flooded gum", "stinking gum", "slaty gum"],

"Juniperus sp.":
    [],

"Ulmus procera":
    ["Englische Ulme", "English elm", "Englische Ulme", "English elm", "English elm", "English elm", "English cork elm", "English elm", "English cork elm", "English elm", "English elm"],

"Tibouchina granulosa":
    ["Brazilian glorytree", "quaresmeira", "Brazilian glorytree"],

"Eucalyptus haemastoma":
    ["scribbly gum", "snappy gum", "white gum", "scribbly gum", "scribbly gum", "scribbly gum", "snappy gum", "white gum", "scribbly gum"],

"Ficus elastica":
    [],

"Albizia julibrissin Durazz":
    ["mimosa", "powderpuff tree", "silk tree", "silktree", ""],

"Morus Sp.":
    [],

"Morus nigra":
    ["Black Mulberry", "black mulberry"],

"Michelia doltsopa":
    ["sweet michelia", "sweet michelia"],

"Thuja sp.":
    [],

"Cupressus macrocarpa":
    ["Monterey Cypress", "Monterey cypress", "Monterey-pine"],

"Pinus sp.":
    [],

"Michelia figo":
    [],

"Cytisus scoparius":
    ["Besenginster", "European broom", "Irish broom", "Scotch broom", "broomtops", "common broom", "genêt à balais", "giesta", "European broom", "Irish broom", "broom", "broomtops", "common broom", "Scotch broom", "Scottish broom", "English broom"],

"Ligustrum ovalifolium Aurea":
    [],

"Eucalyptus longifolia":
    ["woollybutt", "woollybutt", "woollybutt", "woollybutt"],

"Eucalyptus leucoxylon Rosea":
    [],

"Laurus nobilis":
    ["Bay",  "sweet bay", "bay", "bay laurel", "bay-leaf laurel", "Grecian laurel", "laurel",  "Bay Laurel", "Echter Lorbeerbaum",  "sweet bay", "bay laurel", "sweet bay", "yue gui", "gekkeiju", "bay", "bay laurel", "bay-leaf laurel", "Grecian laurel", "louro-comum", "louro-de-apolônio", "louro-europeu", "laurel", "sweet bay"],

"Fraxinus excelsior Aurea":
    [],

"Syzygium jambos":
    ["Malabar plum", "Rosenapfelbaum", "ahi‘a papa‘a", "apel en wai", "fa palangi", "fekika papalangi", "haia", "hehea ha‘amoa", "iouen wai", "iouen wai", "jambos", "jambosier", "jambrosade", "jamrosa", "jamrosa", "jamrosat", "jamrosier", "kavika ni India", "kavika ni vavalangi", "kavika ni vavalangi", "ka‘ika", "ka‘ika papa‘a", "ka‘ika takataka", "ka‘ika varani", "manzana rosa", "pomarrosa", "pomme-rose", "pommier rose", "prunier de Malabar", "rose apple", "rose-apple", "seasea palagi", "yambo", "youenwai", "‘ohi‘a loke"],

"Malus sp.":
    [],

"Eucalyptus cladocalyx":
    ["Sugargum", "sugargum", "sugar gum", "Sugargum", "sugargum", "sugar gum", "sugar gum"],

"Unknown Species":
    [],

"Tibouchina":
    [],

"Ceratonia siliqua":
    ["carob", "locust-bean", "St. John's-bread", "St. John's bread", "carob", "Algaroba", "Carob", "Carob Tree", "Caroubier", "Locust Bean", "Locust Tree", "St. John's Bread", "St.John's Bread", "St. John's bread", "carob", "St. John's-bread", "carob", "carob", "locust-bean", "St. John's-bread"],

"Thuja plicata":
    [],

"Arbutus unedo":
    ["Madrono", "Western Strawberry-tree", "Westlicher Erdbeerbaum", "strawberry tree", "strawberry tree", "strawberry-tree", "Irish strawberry-tree", "arbutus", "Madrono", "Western Strawberry-tree", "Westlicher Erdbeerbaum", "strawberry tree"],

"Acacia longifolia":
    ["Sydney golden wattle", "acácia", "acácia-de-espigas", "acácia-de-folhas-longas", "acácia-marítima", "acácia-trinervis", "golden wattle", "langblaarwattel", "long-leaf wattle", "salgueiro-amarelo", "sallow wattle", "western yarrow", "golden-rods", "longleaf wattle", "sallow wattle", "Sydney golden wattle", "coastal wattle"],

"Caesalpinia ferrea":
    ["Leopard Tree", "Pau-ferro", "Leopard Tree", "Pau-ferro", "Brazilian ironwood", "leopardtree", "Brazilian ironwood", "leopardtree"],

"Persea sp.":
    [],

"Archidendron mullerianum":
    [],

"Ceratopetalum gummiferum":
    ["Christmasbush", "Christmasbush"],

"Callitris columellaris":
    ["Bribie Island Pine", "Cypress Pine", "Murray Pine", "Murray River Pine", "Northern Cypress-pine", "Western Cypress", "Western Sand Cypress", "White Cypress-pine", "White pine", "Bribie Island Pine", "Cypress Pine", "Murray Pine", "Murray River Pine", "Northern Cypress-pine", "Western Cypress", "Western Sand Cypress", "White Cypress-pine", "White pine", "Bribie Island-pine", "Murray River-pine", "Murray-pine", "northern cypress-pine", "sand-cypress", "western-cypress", "cypress-pine", "western cypress-pine", "coast cypress-pine", "coastal-cypress", "western sand-cypress", "white cypress-pine", "white-pine", "slender native cypress-pine", "Bribie Island-pine", "Murray River-pine", "Murray-pine", "northern cypress-pine", "sand-cypress", "western-cypress", "cypress-pine", "western cypress-pine", "coast cypress-pine", "coastal-cypress", "western sand-cypress", "white cypress-pine", "white-pine", "slender native cypress-pine"],

"Eucalyptus camaldulensis x ova":
    [],

"Alphitonia excelsa":
    ["mountain-ash", "pink-almond", "red tweedie", "red-almond", "red-ash", "soaptree", "white-myrtle", "whiteleaf", "Cooper's-wood", "leather-jacket"],

"Zelkova serrata":
    ["Japanese zelkova", "Zelkova tree, Sawleaf zelkova, Japanse zelkova",  ", ", "Japanese zelkova", ";", "Japanese zelkova", "saw-leaf zelkova", "water-elm", "ju shu", "Japanese zelkova", "Japanese-elm", "keyaki", "keyaki", "Japanese zelkova", "saw-leaf zelkova", "water-elm", "Japanese zelkova", "Japanese-elm", "keyaki", ";"],

"Melaleuca bracteata 'Rev. Gold":
    [],

"Eucalyptus leucoxylon":
    ["white ironbark", "yellow gum", "eucalyptus", "white ironbark", "South Australian blue gum", "blue gum", "inland blue gum", "water gum", "yellow gum", "large-fruit blue gum", "large-fruit yellow gum", "red-flower yellow gum", "small-fruit yellow gum", "white ironbark", "white ironbark", "yellow gum", "South Australian blue gum", "blue gum", "inland blue gum", "water gum", "yellow gum", "large-fruit blue gum", "large-fruit yellow gum", "red-flower yellow gum", "small-fruit yellow gum", "white ironbark"],

"Banksia sp.":
    [],

"Magnolia soulangeana":
    ["Chinese magnolia", "Chinese magnolia", "saucer magnolia", "Chinese magnolia", "saucer magnolia", "Chinese magnolia", "Chinese magnolia", "saucer magnolia", "saucer magnolia", "Chinese magnolia"],

"Cedrus deodara":
    ["Deodar cedar", "Himalayan cedar", "himalakski cedar"],

"Banksia ericifolia":
    ["heath-leaf banksia", "heath banksia", "heath-leaf banksia", "heath-leaf banksia", "heath banksia", "heath-leaf banksia"],

"Eucalyptus viminalis":
    ["manna gum", "manna gum", "ribbon gum", "rough-bark manna gum", "white gum", "manna gum", "manna gum", "ribbon gum", "rough-bark manna gum", "white gum"],

"Phoenix sp.":
    [],

"Calodendron capense":
    [],

"Paulownia tomentosa":
    ["empress tree", "foxglove-tree", "karritree", "kiri", "princess tree", "Foxglove tree",   "Blauglockenbaum", "Princesstree", "empress tree", "foxglove-tree", "karritree", "kiri", "princess tree", "foxglovetree", "mao pao tong", "paulovnia-real", "princesstree", "quiri", "karritree", "kiri", "empresstree", "princess tree", "royal paulownia", "foxglovetree", "princesstree", "karritree", "empresstree", "Blauglockenbaum", "Princesstree"],

"Melicope elleryana":
    ["pink doughwood", "pink-euodia", "pink-flower doughwood", "pink-flower-euodia", "pink doughwood", "pink-euodia", "pink-flower doughwood", "pink-flower-euodia"],

"Thevetia peruviana":
    ["Thevetie", "adelfa amarilla", "be still tree", "cabalonga", "chirca", "foreigner's tree", "geel-oleander", "irelepsech", "jacapa", "kanneeta", "koneta", "loandro-amarelo", "luckynut", "nohomalie", "oléandre jaune", "piti", "poupou", "pua", "venevene", "yellow oleander"],

"Davidia involucrata":
    ["handkerchief-tree", "dovetree", "ghost-tree", "handkerchief-tree", "dovetree", "ghost-tree"],

"Pyrus communis":
    [],

"Pinus patula":
    ["patula pine", "Mexican weeping pine", "Mexican yellow pine", "jelecote pine", "Jelecote pine", "Mexican weeping pine", "Spreading-leaved pine", "Mexican weeping pine", "jelecote pine", "pino colorado", "patula pine", "Mexican weeping pine", "Mexican yellow pine", "jelecote pine"],

"Acacia decurrens":
    ["green wattle"],

"Tibouchina urvilleana":
    ["balmane", "doudoul", "glorybush", "griffe du diable", "lasiandra", "lisandra", "pensée malgache", "princess flower", "purple glorytree", "glorybush", "lasiandra", "princess-flower", "princessflower", "balmane", "doudoul", "glorybush", "griffe du diable", "lasiandra", "lisandra", "pensée malgache", "princess flower", "purple glorytree", "glorybush", "lasiandra", "princess-flower", "princessflower", "lasiandra", "princess flower", "purple glory bush", "princess-flower", "purple glorytree", "glorybush", "lasiandra"],

"Lophostemon confertus Variegat":
    [],

"Brachychiton rupestris":
    ["Queensland bottletree", "narrow-leaf bottletree", "Queensland bottletree", "Queensland rattletree", "narrow-leaf bottletree", "bottletree", "Queensland bottletree", "Queensland rattletree", "narrow-leaf bottletree", "bottletree"],

"Casuarina torulosa":
    [],

"Rothmannia globosa":
    [],

"Quercus robur":
    ["pedunculate oak", "English oak", "English oak", "pedunculate oak", "truffle oak", "dub &ccaron;ereš&ccaron;atyj", "pedunculate oak", "English oak", "xia li", "European oak", "Pedunculate oak", "English oak", "European oak", "pedunculate oak"],

"Eucalyptus linearis":
    [],

"Erythrina crista-galli":
    ["Cockspur Coral Tree", "Coral", "Coral Tree", "Seibo", ";", "crybabytree", "ceibo", "cockspur coraltree", "cry-baby tree", "cockspur coraltree", "cry-baby-tree", "cockspur coraltree", "cry-baby-tree", "crybabytree"],
}

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
syd_map = syd.plot(color="ghostwhite", edgecolor="black")
plt.xlim((151.17,151.28))
plt.ylim((-33.925, -33.845))
# print(plt.axis())
for i, x in enumerate(tree_counts.iteritems()):
    tree = x[0]
    count = x[1]
    temp_df = gdf[gdf.species == tree]
    try:
        common_name = "—"+tree_common_names[tree][0].title()
    except:
        common_name = ""
    label = r"$\it{" + tree + "}$" + f"{common_name} ({count})" if count >= 50 else None
    temp_df.plot(
        ax=syd_map, 
        color=idents[i]["c"], 
        marker=idents[i]["m"], 
        alpha=0.4, 
        label=label, 
        markersize=3
    )

plt.legend(
    loc="upper right", 
    markerscale=2, 
    prop={'size': 6},
    labelspacing=0
)
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
    "all_tree_map", bbox_inches="tight", dpi=600, metadata=metadata, format="png"
)
plt.savefig(
    "all_tree_map", bbox_inches="tight", dpi=600, metadata=metadata, format="svg"
)

#%%
# import contextily as ctx
# tdf = geopandas.read_file(geopandas.datasets.get_path('nybb'))
# ax = tdf.plot(figsize=(10, 10), alpha=0.5, edgecolor='k')
# ctx.add_basemap(ax, url=ctx.providers.Stamen.TonerLite)
# ax.set_axis_off()
# %%
