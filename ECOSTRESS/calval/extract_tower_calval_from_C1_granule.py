from logging import getLogger
import warnings

import numpy as np
import pandas as pd
import rasters as rt

import colored_logging as cl

from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule

from ECOSTRESS.calval.tower import Tower
from ECOSTRESS.calval.generate_3x3_subset_geometry import generate_3x3_subset_geometry

logger = getLogger(__name__)

def extract_tower_calval_from_C1_granule(
        C1_granule: ECOSTRESSSwathGranule, 
        tower: Tower,
        subset_3x3: rt.RasterGrid = None) -> pd.DataFrame:
    calval_rows = []

    tower_name = tower.name
    lat = tower.lat
    lon = tower.lon
    tower_point_latlon = rt.Point(lon, lat)

    if subset_3x3 is None:
        subset_3x3 = generate_3x3_subset_geometry(C1_granule.geometry, tower_point_latlon)

    for variable in sorted(set(C1_granule.variables)):
        if "quality" in variable:
            continue

        try:
            variable_subset = C1_granule.variable(variable, geometry=subset_3x3)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                median_value = np.nanmedian(variable_subset)

            if np.isnan(median_value):
                continue

            orbit = C1_granule.orbit
            scene = C1_granule.scene
            tile = ""

            FieldOfViewObstruction = ""
            OrbitCorrectionPerformed = ""
            QAPercentCloudCover = ""
            QAPercentGoodQuality = ""
            filename = C1_granule.filename

            logger.info(f"{cl.name(tower)} {cl.time(C1_granule.time_UTC)} {cl.name(variable)}: {median_value}")

            calval_rows.append([
                tower_name,
                orbit,
                scene,
                tile,
                C1_granule.time_UTC,
                C1_granule.product,
                variable,
                median_value,
                1,
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
        "FieldOfViewObstruction",
        "OrbitCorrectionPerformed",
        "QAPercentCloudCover",
        "QAPercentGoodQuality",
        "filename",
        "lat",
        "lon"
    ])

    return df
