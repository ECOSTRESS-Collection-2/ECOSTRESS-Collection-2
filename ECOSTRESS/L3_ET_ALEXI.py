"""
This module handles ingest of an ECOSTRESS L3 ET PT-JPL granule.
Developed by Gregory Halverson at the Jet Propulsion Laboratory
"""
import json
import logging
from os.path import exists
from typing import List

from ECOSTRESS_colors import ET_COLORMAP

from rasters import Raster

from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule
from ECOSTRESS.ECOSTRESS_gridded_tiled_granule import ECOSTRESSGriddedGranule, ECOSTRESSTiledGranule

__author__ = "Gregory Halverson"


class L3ETALEXI(ECOSTRESSSwathGranule):
    logger = logging.getLogger(__name__)

    _DATA_GROUP = "EVAPOTRANSPIRATION ALEXI"
    _PRODUCT_METADATA_GROUP = "L3_ET_ALEXI Metadata"

    def __init__(
            self,
            L3_ET_ALEXI_filename: str,
            L2_CLOUD_filename: str = None,
            L1B_GEO_filename: str = None,
            **kwargs):
        ECOSTRESSSwathGranule.__init__(
            self,
            product_filename=L3_ET_ALEXI_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            data_group_name=self._DATA_GROUP,
            product_metadata_group_name=self._PRODUCT_METADATA_GROUP,
            **kwargs
        )

        if not exists(L3_ET_ALEXI_filename):
            raise IOError(f"L3 ALEXI file does not exist: {L3_ET_ALEXI_filename}")
    @property
    def ETdaily(self) -> Raster:
        return self.variable("ETdaily")

    def __repr__(self):
        display_dict = {
            "L3 disALEXI": self.product_filename,
            "L2 CLOUD:": self.L2_CLOUD_filename,
            "L1B GEO:": self.L1B_GEO_filename
        }

        display_string = f"L3 disALEXI Granule\n{json.dumps(display_dict, indent=2)}"

        return display_string


class L3TETALEXI(ECOSTRESSTiledGranule, L3ETALEXI):
    _PRIMARY_VARIABLE = "ETdaily"
    _PRODUCT_NAME = "L3T_ET_ALEXI"

    @property
    def primary_variable(self) -> str:
        return "ETdaily"

    @property
    def variables(self) -> List[str]:
        return ["ETdaily", "ETdailyUncertainty"]

class L3GETALEXI(ECOSTRESSGriddedGranule, L3ETALEXI):
    _TILE_CLASS = L3TETALEXI
    _GRID_NAME = "ECO_L3G_ET_ALEXI_70m"
    _PRIMARY_VARIABLE = "ETdaily"
    _GRANULE_PREVIEW_CMAP = ET_COLORMAP

    VARIABLE_NAMES = [
        "ETdaily",
        "ETdailyUncertainty"
    ]

    def __init__(self, L3G_ET_ALEXI_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L3G_ET_ALEXI_filename)
