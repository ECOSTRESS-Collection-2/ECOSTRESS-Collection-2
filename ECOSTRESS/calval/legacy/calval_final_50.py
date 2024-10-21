import logging
import warnings
from glob import glob
from os.path import join
import pandas as pd
import rasters as rt
import numpy as np
import cl
from ECOSTRESS import open_granule

logger = logging.getLogger(__name__)

# product_directory = "/project/sandbox/ECOv002_L3T_L4T/output"
product_directory = "/project/sandbox/calval_v4/output"

tower_locations = pd.read_csv("final_50_sites.csv")

pattern = join(product_directory, "**", "*.zip")
logging.info(f"searching product files: {cl.val(pattern)}")
filenames = sorted(glob(pattern, recursive=True))
logger.info(f"found {cl.val(len(filenames))} files")

calval_rows = []

for filename in filenames:
    logger.info(f"file: {cl.file(filename)}")

    try:
        granule = open_granule(filename)
    except Exception as e:
        logger.exception(e)
        logger.warning(f"unable to load file: {filename}")
        continue

    for i, tower_locations_row in tower_locations.iterrows():
        tower = tower_locations_row.tower
        lat = tower_locations_row.latitude
        lon = tower_locations_row.longitude
        tower_point_latlon = rt.Point(lon, lat)
        # vegetation = tower_locations_row.vegetation
        # climate = tower_locations_row.climate
        # elevation = tower_locations_row.elevation
        # MAT_C = tower_locations_row.MAT_C
        # MAP_mm = tower_locations_row.MAP_mm

        if granule.geometry.intersects(tower_point_latlon.to_crs(granule.geometry.crs)):
            tower_row, tower_col = granule.geometry.index_point(tower_point_latlon.to_crs(granule.geometry.crs))
            rows, cols = granule.geometry.shape
            subset_3x3 = granule.geometry[max(tower_row - 1, 0):min(tower_row + 2, rows - 1),
                         max(tower_col - 1, 0):min(tower_col + 2, cols - 1)]

            for variable in sorted(set(granule.variables) - {"water", "cloud", "QC"}):
                if "quality" in variable:
                    continue

                try:
                    variable_subset = granule.variable(variable, geometry=subset_3x3)

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        median_value = np.nanmedian(variable_subset)

                    if np.isnan(median_value):
                        continue

                    orbit = granule.orbit
                    scene = granule.scene
                    tile = granule.tile

                    logger.info(f"{cl.name(tower)} {cl.time(granule.time_UTC)} {cl.name(variable)}: {median_value}")

                    calval_rows.append([
                        tower,
                        orbit,
                        scene,
                        tile,
                        granule.time_UTC,
                        granule.product,
                        variable,
                        median_value,
                        # vegetation,
                        # climate,
                        # MAT_C,
                        # MAP_mm,
                        # elevation,
                        lat,
                        lon
                    ])
                except Exception as e:
                    logger.exception(e)
                    logger.warning(f"unable to read {cl.name(tower)} {cl.time(granule.time_UTC)} {cl.name(variable)}")
                    continue

df = pd.DataFrame(calval_rows, columns=[
    "tower",
    "orbit",
    "scene",
    "tile",
    "time_UTC",
    "product",
    "variable",
    "med3x3",
    # "vegetation",
    # "climate",
    # "MAT_C",
    # "MAP_mm",
    # "elevation",
    "lat",
    "lon"
])

df.to_csv("/project/sandbox/ECOv002_L3T_L4T/calval_final_50.csv")
