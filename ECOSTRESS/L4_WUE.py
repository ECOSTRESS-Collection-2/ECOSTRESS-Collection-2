"""
This module handles ingest of an ECOSTRESS L4 WUE granule.
Developed by Gregory Halverson at the Jet Propulsion Laboratory
"""
import json
import logging
from os.path import exists

from rasters import Raster

from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule
from ECOSTRESS.ECOSTRESS_gridded_tiled_granule import ECOSTRESSGriddedGranule, ECOSTRESSTiledGranule

__author__ = "Gregory Halverson"

from ECOSTRESS_colors import NDVI_COLORMAP, GPP_COLORMAP

GRANULE_PREVIEW_CMAP = NDVI_COLORMAP

class L4WUE(ECOSTRESSSwathGranule):
    logger = logging.getLogger(__name__)

    _DATA_GROUP = "Water Use Efficiency"
    _PRODUCT_METADATA_GROUP = "L4_WUE_Metadata"
    _GRANULE_PREVIEW_CMAP = GRANULE_PREVIEW_CMAP

    def __init__(
            self,
            L4_WUE_filename: str,
            L2_CLOUD_filename: str = None,
            L1B_GEO_filename: str = None,
            **kwargs):
        super().__init__(
            product_filename=L4_WUE_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            data_group_name=self._DATA_GROUP,
            product_metadata_group_name=self._PRODUCT_METADATA_GROUP,
            **kwargs
        )

        if not exists(L4_WUE_filename):
            raise IOError(f"L4 WUE file does not exist: {L4_WUE_filename}")

    def __repr__(self):
        display_dict = {
            "L4 WUE": self.product_filename,
            "L2 CLOUD:": self.L2_CLOUD_filename,
            "L1B GEO:": self.L1B_GEO_filename
        }

        display_string = f"L4 WUE Granule\n{json.dumps(display_dict, indent=2)}"

        return display_string

class ECOv001L4WUEPTJPL(L4WUE):
    @property
    def WUE(self) -> Raster:
        return self.variable("WUEavg", cmap=GPP_COLORMAP)

class ECOv002L4WUE(L4WUE):
    @property
    def WUE(self) -> Raster:
        return self.variable("WUE", cmap=GPP_COLORMAP)

    @property
    def GPP(self) -> Raster:
        return self.variable("GPP", cmap=GPP_COLORMAP)

    @property
    def GPP_umol(self) -> Raster:
        return self.GPP

    @property
    def GPP_mol(self) -> Raster:
        return self.GPP_umol / 1000000

    @property
    def GPP_g(self) -> Raster:
        GPP_g = self.GPP_mol * 12.011

        return GPP_g

class L4TWUE(ECOSTRESSTiledGranule, ECOv002L4WUE):
    _PRIMARY_VARIABLE = "WUE"
    _PRODUCT_NAME = "L4T_WUE"
    _GRANULE_PREVIEW_CMAP = GPP_COLORMAP


class L4GWUE(ECOSTRESSGriddedGranule, ECOv002L4WUE):
    _TILE_CLASS = L4TWUE
    _GRID_NAME = "ECO_L4G_WUE_70m"
    _GRANULE_PREVIEW_CMAP = GRANULE_PREVIEW_CMAP
    _PRIMARY_VARIABLE = "WUE"

    _DATASET_NAME_TRANSLATIONS = {
        "WUEavg": "WUE"
    }

    VARIABLE_NAMES = [
        "WUE",
        "GPP"
    ]

    def __init__(self, L4G_WUE_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L4G_WUE_filename)
