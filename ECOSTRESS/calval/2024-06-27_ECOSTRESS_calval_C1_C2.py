import logging
import warnings
from glob import glob
from os.path import join, abspath, dirname
import pandas as pd
import rasters as rt
import numpy as np
import cl
from ECOSTRESS.find_ECOSTRESS_C1_scene import find_ECOSTRESS_C1_scene
from ECOSTRESS import open_granule

logger = logging.getLogger(__name__)

C2_product_directory = "/project/sandbox/ECOv002_L3T_L4T/output"
# C2_product_directory = "/project/sandbox/calval_v4/output"

tower_locations = pd.read_csv(join(abspath(dirname(__file__)), "2024-06-27_ECOSTRESS_calval_sites.csv"))

C2_search_pattern = join(C2_product_directory, "**", "*.zip")
logging.info(f"searching product files: {cl.val(C2_search_pattern)}")
C2_filenames = sorted(glob(C2_search_pattern, recursive=True))

logger.info(f"found {cl.val(len(C2_filenames))} files")

C1_PRODUCT_NAMES = {
    "L3T_JET": "L3_ET_PT-JPL",
    "L3T_ET_ALEXI": "L3_ET_ALEXI",
    "L4T_ESI": "L4_ESI_PT-JPL",
    "L4T_ESI_ALEXI": "L4_ESI_ALEXI",
    "L4T_WUE": "L4_WUE"
}

calval_rows = []

for C2_filename in C2_filenames:
    logger.info(f"file: {cl.file(C2_filename)}")

    try:
        C2_granule = open_granule(C2_filename)
    except Exception as e:
        logger.exception(e)
        logger.warning(f"unable to load file: {C2_filename}")
        continue

    for i, tower_locations_row in tower_locations.iterrows():
        tower = tower_locations_row.tower
        lat = tower_locations_row.latitude
        lon = tower_locations_row.longitude
        tower_point_latlon = rt.Point(lon, lat)


        if C2_granule.geometry.intersects(tower_point_latlon.to_crs(C2_granule.geometry.crs)):
            tower_row, tower_col = C2_granule.geometry.index_point(tower_point_latlon.to_crs(C2_granule.geometry.crs))
            rows, cols = C2_granule.geometry.shape
            subset_3x3 = C2_granule.geometry[max(tower_row - 1, 0):min(tower_row + 2, rows - 1),
                         max(tower_col - 1, 0):min(tower_col + 2, cols - 1)]

            for variable in sorted(set(C2_granule.variables) - {"water", "cloud", "QC"}):
                if "quality" in variable:
                    continue

                try:
                    variable_subset = C2_granule.variable(variable, geometry=subset_3x3)

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        median_value = np.nanmedian(variable_subset)

                    if np.isnan(median_value):
                        continue

                    orbit = C2_granule.orbit
                    scene = C2_granule.scene
                    tile = C2_granule.tile

                    logger.info(f"{cl.name(tower)} {cl.time(C2_granule.time_UTC)} {cl.name(variable)}: {median_value}")

                    calval_rows.append([
                        tower,
                        orbit,
                        scene,
                        tile,
                        C2_granule.time_UTC,
                        C2_granule.product,
                        variable,
                        median_value,
                        2,
                        lat,
                        lon
                    ])
                except Exception as e:
                    logger.exception(e)
                    logger.warning(f"unable to read {cl.name(tower)} {cl.time(C2_granule.time_UTC)} {cl.name(variable)}")
                    continue
    
            if C2_granule.product in C1_PRODUCT_NAMES:
                C1_product = C1_PRODUCT_NAMES[C2_granule.product]
                
                C1_filename_dict = find_ECOSTRESS_C1_scene(
                    orbit=C2_granule.orbit,
                    scene=C2_granule.scene,
                    date_UTC=C2_granule.date_UTC
                )

                if C1_product in C1_filename_dict:
                    C1_filename = C1_filename_dict[C1_product]
                    C1_L1B_GEO_filename = C1_filename_dict["L1B_GEO"]
                    C1_L2_cloud_filename = C1_filename_dict["L2_CLOUD"]
                    
                    C1_granule = open_granule(
                        filename=C1_filename,
                        L1B_GEO_filename=C1_L1B_GEO_filename,
                        L2_CLOUD_filename=C1_L2_cloud_filename
                    )

                    tower_row, tower_col = C2_granule.geometry.index_point(tower_point_latlon.to_crs(C2_granule.geometry.crs))
                    rows, cols = C2_granule.geometry.shape
                    row_min = max(tower_row - 1, 0)
                    row_max = min(tower_row + 2, rows - 1)
                    col_min = max(tower_col - 1, 0)
                    col_max = min(tower_col + 2, cols - 1)
                    subset_3x3 = C2_granule.geometry[row_min:row_max, col_min:col_max]

                    for variable in sorted(set(C1_granule.variables)):
                        if "quality" in variable:
                            continue

                        try:
                            # variable_subset = C1_granule.variable(variable, geometry=subset_3x3)
                            variable_full = C1_granule.variable(variable)
                            variable_subset = variable_full[row_min:row_max, col_min:col_max]

                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore")
                                median_value = np.nanmedian(variable_subset)

                            if np.isnan(median_value):
                                continue

                            orbit = C1_granule.orbit
                            scene = C1_granule.scene
                            tile = ""

                            logger.info(f"{cl.name(tower)} {cl.time(C1_granule.time_UTC)} {cl.name(variable)}: {median_value}")

                            calval_rows.append([
                                tower,
                                orbit,
                                scene,
                                tile,
                                C1_granule.time_UTC,
                                C1_granule.product,
                                variable,
                                median_value,
                                1,
                                lat,
                                lon
                            ])
                        except Exception as e:
                            logger.exception(e)
                            logger.warning(f"unable to read {cl.name(tower)} {cl.time(C1_granule.time_UTC)} {cl.name(variable)}")
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
    "collection",
    "lat",
    "lon"
])

df.to_csv("/project/sandbox/ECOv002_L3T_L4T/2024-06-27_ECOSTRESS_calval_C1_C2.csv")
