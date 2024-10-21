from logging import getLogger
import warnings

import numpy as np
import pandas as pd
import rasters as rt

import colored_logging as cl

from ECOSTRESS.ECOSTRESS_gridded_tiled_granule import ECOSTRESSTiledGranule

from ECOSTRESS.calval.tower import Tower
from ECOSTRESS.calval.generate_3x3_subset_geometry import generate_3x3_subset_geometry

logger = getLogger(__name__)

def extract_tower_calval_from_C2_granule(
        C2_granule: ECOSTRESSTiledGranule, 
        tower: Tower,
        subset_3x3: rt.RasterGrid = None) -> pd.DataFrame:
    calval_rows = []

    tower_name = tower.name
    lat = tower.lat
    lon = tower.lon
    tower_point_latlon = rt.Point(lon, lat)

    if subset_3x3 is None:
        subset_3x3 = generate_3x3_subset_geometry(C2_granule.geometry, tower_point_latlon)

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

            if "FieldOfViewObstruction" in C2_granule.standard_metadata:
                FieldOfViewObstruction = C2_granule.standard_metadata["FieldOfViewObstruction"]
            else:
                FieldOfViewObstruction = ""
            
            if "OrbitCorrectionPerformed" in C2_granule.product_metadata:
                OrbitCorrectionPerformed = C2_granule.product_metadata["OrbitCorrectionPerformed"]
            else:
                OrbitCorrectionPerformed = ""

            if "QAPercentCloudCover" in C2_granule.product_metadata:
                QAPercentCloudCover = C2_granule.product_metadata["QAPercentCloudCover"]
            else:
                QAPercentCloudCover = ""

            if "QAPercentGoodQuality" in C2_granule.product_metadata:
                QAPercentGoodQuality = C2_granule.product_metadata["QAPercentGoodQuality"]
            else:
                QAPercentGoodQuality = ""

            filename = C2_granule.filename

            logger.info(f"{cl.name(tower)} {cl.time(C2_granule.time_UTC)} {cl.name(variable)}: {median_value}")

            calval_rows.append([
                tower_name,
                orbit,
                scene,
                tile,
                C2_granule.time_UTC,
                C2_granule.product,
                variable,
                median_value,
                2,
                FieldOfViewObstruction,
                OrbitCorrectionPerformed,
                QAPercentCloudCover,
                QAPercentGoodQuality,
                filename,
                lat,
                lon
            ])
        except Exception as e:
            logger.exception(e)
            logger.warning(f"unable to read {cl.name(tower)} {cl.time(C2_granule.time_UTC)} {cl.name(variable)}")
            continue

    calval_rows

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
        "FieldOfViewObstruction",
        "OrbitCorrectionPerformed",
        "QAPercentCloudCover",
        "QAPercentGoodQuality",
        "filename",
        "lat",
        "lon"
    ])

    return df
