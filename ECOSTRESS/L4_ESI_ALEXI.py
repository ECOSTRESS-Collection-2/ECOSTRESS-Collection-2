"""
This module handles ingest of an ECOSTRESS L4 ESI disALEXI granule.
Developed by Gregory Halverson at the JESI Propulsion Laboratory
"""
import json
import logging
from os.path import exists
from typing import List

from ECOSTRESS_colors import ET_COLORMAP

from ECOSTRESS.ECOSTRESS_gridded_tiled_granule import ECOSTRESSGriddedGranule, ECOSTRESSTiledGranule
from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule

__author__ = "Gregory Halverson"

from rasters import RasterGeometry

PRIMARY_VARIABLE = "ESIdaily"

logger = logging.getLogger(__name__)


class L4ESIALEXI(ECOSTRESSSwathGranule):
    _DATA_GROUP = "EVAPORATIVE STRESS INDEX ALEXI"
    _PRODUCT_METADATA_GROUP = "L4_ESI_ALEXI Metadata"
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE

    def __init__(
            self,
            L4_ESI_ALEXI_filename: str,
            L2_CLOUD_filename: str = None,
            L1B_GEO_filename: str = None,
            **kwargs):
        super().__init__(
            product_filename=L4_ESI_ALEXI_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            data_group_name=self._DATA_GROUP,
            product_metadata_group_name=self._PRODUCT_METADATA_GROUP,
            **kwargs
        )

        if not exists(L4_ESI_ALEXI_filename):
            raise IOError(f"L4 ESI disALEXI file does not exist: {L4_ESI_ALEXI_filename}")

    def __repr__(self):
        display_dict = {
            "L4 ESI disALEXI": self.product_filename,
            "L2 CLOUD:": self.L2_CLOUD_filename,
            "L1B GEO:": self.L1B_GEO_filename
        }

        display_string = f"L4 ESI disALEXI Granule\n{json.dumps(display_dict, indent=2)}"

        return display_string


class L4TESIALEXI(ECOSTRESSTiledGranule):
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE
    _PRODUCT_NAME = "L4T_ESI_ALEXI"

    # @property
    # def geometry(self) -> RasterGeometry:
    #     if self._geometry is not None:
    #         return self._geometry
    #     else:
    #         geometry = self.variable("ESIdaily").geometry
    #         self._geometry = geometry
    #
    #         return geometry

    @property
    def primary_variable(self) -> str:
        return "ESIdaily"

    @property
    def variables(self) -> List[str]:
        return ["ESIdaily", "ESIdailyUncertainty"]


class L4GESIALEXI(ECOSTRESSGriddedGranule, L4ESIALEXI):
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE
    _TILE_CLASS = L4TESIALEXI
    _GRID_NAME = "ECO_L4G_ESI_ALEXI_70m"
    _GRANULE_PREVIEW_CMAP = ET_COLORMAP

    VARIABLE_NAMES = [
        "ESIdaily",
        "ESIdailyUncertainty"
    ]

    def __init__(self, L4G_ESI_ALEXI_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L4G_ESI_ALEXI_filename)
