"""
This module handles ingest of an ECOSTRESS L4 ESI granule.
Developed by Gregory Halverson at the JESI Propulsion Laboratory
"""
import json
import logging
from os.path import exists

from ECOSTRESS_colors import ET_COLORMAP
from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule
from ECOSTRESS.ECOSTRESS_gridded_tiled_granule import ECOSTRESSGriddedGranule, ECOSTRESSTiledGranule

__author__ = "Gregory Halverson"

GRANULE_PREVIEW_CMAP = ET_COLORMAP

class L4ESI(ECOSTRESSSwathGranule):
    logger = logging.getLogger(__name__)

    _DATA_GROUP = "Evaporative Stress Index PT-JPL"
    _PRODUCT_METADATA_GROUP = "L4_ESI Metadata"
    _GRANULE_PREVIEW_CMAP = GRANULE_PREVIEW_CMAP

    def __init__(
            self,
            L4_ESI_filename: str,
            L2_CLOUD_filename: str = None,
            L1B_GEO_filename: str = None,
            **kwargs):
        super().__init__(
            product_filename=L4_ESI_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            data_group_name=self._DATA_GROUP,
            product_metadata_group_name=self._PRODUCT_METADATA_GROUP,
            **kwargs
        )

        if not exists(L4_ESI_filename):
            raise IOError(f"L4 ESI file does not exist: {L4_ESI_filename}")

    def __repr__(self):
        display_dict = {
            "L4 ESI": self.product_filename,
            "L2 CLOUD:": self.L2_CLOUD_filename,
            "L1B GEO:": self.L1B_GEO_filename
        }

        display_string = f"L4 ESI Granule\n{json.dumps(display_dict, indent=2)}"

        return display_string


class L4TESI(ECOSTRESSTiledGranule):
    _PRIMARY_VARIABLE = "ESI"
    _PRODUCT_NAME = "L4T_ESI"


class L4GESI(ECOSTRESSGriddedGranule, L4ESI):
    _TILE_CLASS = L4TESI
    _GRID_NAME = "ECO_L4G_ESI_70m"
    _GRANULE_PREVIEW_CMAP = GRANULE_PREVIEW_CMAP
    _PRIMARY_VARIABLE = "ESI"

    VARIABLE_NAMES = [
        "ESI",
        "PET"
    ]

    _DATASET_NAME_TRANSLATIONS = {
        "ESIavg": "ESI"
    }

    def __init__(self, L4G_ESI_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L4G_ESI_filename)
