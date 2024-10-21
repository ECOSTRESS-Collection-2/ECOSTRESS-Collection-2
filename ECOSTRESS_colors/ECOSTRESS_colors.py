from matplotlib.colors import LinearSegmentedColormap

DEFAULT_COLORMAP = "jet"

ET_COLORMAP = LinearSegmentedColormap.from_list("ET", [
    "#f6e8c3",
    "#d8b365",
    "#99974a",
    "#53792d",
    "#6bdfd2",
    "#1839c5"
])

SM_COLORMAP = LinearSegmentedColormap.from_list("SM", [
    "#f6e8c3",
    "#d8b365",
    "#99894a",
    "#2d6779",
    "#6bdfd2",
    "#1839c5"
])

RH_COLORMAP = SM_COLORMAP

NDVI_COLORMAP_ABSOLUTE = LinearSegmentedColormap.from_list(
    name="NDVI",
    colors=[
        (0, "#0000ff"),
        (0.4, "#000000"),
        (0.5, "#745d1a"),
        (0.6, "#e1dea2"),
        (0.8, "#45ff01"),
        (1, "#325e32")
    ]
)

# NDVI_COLORMAP = LinearSegmentedColormap.from_list(
#     name="NDVI",
#     colors=[
#         "#0000ff",
#         "#000000",
#         "#745d1a",
#         "#e1dea2",
#         "#45ff01",
#         "#325e32"
#     ]
# )

NDVI_COLORMAP = NDVI_COLORMAP_ABSOLUTE

ALBEDO_COLORMAP = LinearSegmentedColormap.from_list(name="albedo", colors=["black", "white"])
CLOUD_CMAP = LinearSegmentedColormap.from_list(name="cloud", colors=["black", "white"])

RN_COLORMAP = "jet"
TA_COLORMAP = "jet"

WATER_COLORMAP = LinearSegmentedColormap.from_list(name="water", colors=["black", "#0eeded"])
CLOUD_COLORMAP = LinearSegmentedColormap.from_list(name="water", colors=["black", "white"])

# GPP_COLORMAP = LinearSegmentedColormap.from_list(
#     name="GPP",
#     colors=[
#         "#000000",
#         "#325e32"
#     ]
# )

# GPP_COLORMAP = LinearSegmentedColormap.from_list(
#     name="GPP",
#     colors=[
#         "#000000",
#         "#636821",
#         "#325e32",
#         "#a6ff01",
#         "#00ff00"
#     ]
# )

GPP_COLORMAP = LinearSegmentedColormap.from_list(
    name="GPP",
    colors=[
        "#000000",
        "#bdae08",
        "#325e32",
        "#a6ff01",
        "#00ff00"
    ]
)
