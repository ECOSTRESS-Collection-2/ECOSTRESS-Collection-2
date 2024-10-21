import logging
import warnings
from glob import glob
from os.path import join, abspath, dirname
import pandas as pd
import rasters as rt
import numpy as np
import colored_logging as cl

from ECOSTRESS.find_ECOSTRESS_C1_scene import find_ECOSTRESS_C1_scene
from ECOSTRESS import open_granule

from ECOSTRESS.calval.tower import Tower
from ECOSTRESS.calval.extract_tower_calval_from_C2_granule import extract_tower_calval_from_C2_granule
from ECOSTRESS.calval.extract_tower_calval_from_C1_granule import extract_tower_calval_from_C1_granule

logger = logging.getLogger(__name__)

C2_product_directory = "/project/sandbox/ECOv002_L3T_L4T/output"

tower_locations = pd.read_csv(join(abspath(dirname(__file__)), "2024-08-02_tower_short_list.csv"))

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

calval_df = pd.DataFrame([], columns=[
        "tower",
        "orbit",
        "scene",
        "tile",
        "time_UTC",
        "product",
        "variable",
        "med3x3",
        "collection",
        "FieldOfViewObstruction",
        "OrbitCorrectionPerformed",
        "QAPercentCloudCover",
        "QAPercentGoodQuality",
        "filename",
        "lat",
        "lon"
    ])

for C2_filename in C2_filenames:
    logger.info(f"file: {cl.file(C2_filename)}")

    try:
        C2_granule = open_granule(C2_filename)
    except Exception as e:
        logger.exception(e)
        logger.warning(f"unable to load file: {C2_filename}")
        continue

    tile = C2_granule.tile
    tower_locations_at_tile = tower_locations[tower_locations.tile == tile]

    for i, tower_locations_row in tower_locations_at_tile.iterrows():
        tower_name = tower_locations_row.tower
        lat = tower_locations_row.latitude
        lon = tower_locations_row.longitude
        tower_point_latlon = rt.Point(lon, lat)

        tower = Tower(tower_name, lat, lon)

        if C2_granule.geometry.intersects(tower_point_latlon.to_crs(C2_granule.geometry.crs)):
            C2_df = extract_tower_calval_from_C2_granule(C2_granule, tower)
            calval_df = pd.concat([calval_df, C2_df])
    
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

                    C1_df = extract_tower_calval_from_C1_granule(C1_granule, tower)
                    calval_df = pd.concat([calval_df, C1_df])

calval_df.to_csv("/project/sandbox/ECOv002_L3T_L4T/2024-08-16_ECOSTRESS_calval_C1_C2.csv")
