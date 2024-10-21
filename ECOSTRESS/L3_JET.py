"""
This module handles ingest of an ECOSTRESS L3 ET PT-JPL granule.
Developed by Gregory Halverson at the Jet Propulsion Laboratory
"""
import json
import logging
from os.path import exists

from ECOSTRESS_colors import ET_COLORMAP, SM_COLORMAP, RN_COLORMAP, TA_COLORMAP
from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule
from ECOSTRESS.ECOSTRESS_gridded_tiled_granule import ECOSTRESSGriddedGranule, ECOSTRESSTiledGranule
from ECOSTRESS.ECOSTRESS_granule import ECOSTRESSGranule

__author__ = "Gregory Halverson"

from rasters import Raster

GRANULE_PREVIEW_CMAP = ET_COLORMAP

class L3JET(ECOSTRESSSwathGranule):
    logger = logging.getLogger(__name__)

    _DATA_GROUP_NAME = "EVAPOTRANSPIRATION PT-JPL"
    _PRODUCT_METADATA_GROUP = "L3_ET_PT-JPL Metadata"
    _GRANULE_PREVIEW_CMAP = GRANULE_PREVIEW_CMAP

    def __init__(
            self,
            L3_JET_filename: str,
            L2_CLOUD_filename: str = None,
            L1B_GEO_filename: str = None,
            **kwargs):
        super().__init__(
            product_filename=L3_JET_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            data_group_name=self._DATA_GROUP_NAME,
            product_metadata_group_name=self._PRODUCT_METADATA_GROUP,
            **kwargs
        )

        if not exists(L3_JET_filename):
            raise IOError(f"L3 PTJPLSM file does not exist: {L3_JET_filename}")

    def __repr__(self):
        display_dict = {
            "L3 PT-JPL": self.product_filename,
            "L2 CLOUD:": self.L2_CLOUD_filename,
            "L1B GEO:": self.L1B_GEO_filename
        }

        display_string = f"L3 PT-JPL Granule\n{json.dumps(display_dict, indent=2)}"


        return display_string

    @property
    def STICinst(self) -> Raster:
        return self.variable("STICinst", cmap=ET_COLORMAP)

    @property
    def STICinst_kg(self) -> Raster:
        return self.STICinst / 2450000.0

    @property
    def PTJPLSMinst(self) -> Raster:
        return self.variable("PTJPLSMinst", cmap=ET_COLORMAP)

    @property
    def PTJPLSMinst_kg(self) -> Raster:
        return self.PTJPLSMinst / 2450000.0

    @property
    def MOD16inst(self) -> Raster:
        return self.variable("MOD16inst", cmap=ET_COLORMAP)

    @property
    def MOD16inst_kg(self) -> Raster:
        return self.MOD16inst / 2450000.0

    @property
    def BESSinst(self) -> Raster:
        return self.variable("BESSinst", cmap=ET_COLORMAP)

    @property
    def BESSinst_kg(self) -> Raster:
        return self.BESSinst / 2450000.0

class L3MET(ECOSTRESSGranule):
    pass

class L3SEB(ECOSTRESSGranule):
    pass

class L3SM(ECOSTRESSGranule):
    pass

class L3TBESS(ECOSTRESSTiledGranule):
    _PRIMARY_VARIABLE = "GPP"
    _PRODUCT_NAME = "L3T_BESS"

class L3TJET(ECOSTRESSTiledGranule, L3JET):
    _PRIMARY_VARIABLE = "PTJPLSMinst"
    _PRODUCT_NAME = "L3T_JET"

class L3TMET(ECOSTRESSTiledGranule):
    _PRIMARY_VARIABLE = "Ta"
    _PRODUCT_NAME = "L3T_MET"

class L3TSEB(ECOSTRESSTiledGranule):
    _PRIMARY_VARIABLE = "Rn"
    _PRODUCT_NAME = "L3T_SEB"

class L3TSM(ECOSTRESSTiledGranule):
    _PRIMARY_VARIABLE = "SM"
    _PRODUCT_NAME = "L3T_SM"

    @property
    def SM(self) -> Raster:
        return self.variable("SM")

class L3GJET(ECOSTRESSGriddedGranule, L3JET):
    _TILE_CLASS = L3TJET
    _GRID_NAME = "ECO_L3G_JET_70m"
    _GRANULE_PREVIEW_CMAP = GRANULE_PREVIEW_CMAP
    _PRIMARY_VARIABLE = "PTJPLSMinst"

    VARIABLE_NAMES = [
        "PTJPLSMinst",
        "STICinst",
        "STICcanopy",
        "BESSinst",
        "MOD16inst",
        "ETinstUncertainty",
        "ETdaily",
        "PTJPLSMcanopy",
        "PTJPLSMsoil",
        "PTJPLSMinterception",
        "water",
        "cloud"
    ]

    def __init__(self, L3G_JET_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L3G_JET_filename)

class L3GMET(ECOSTRESSGriddedGranule, L3MET):
    _TILE_CLASS = L3TMET
    _GRID_NAME = "ECO_L3G_MET_70m"
    _GRANULE_PREVIEW_CMAP = TA_COLORMAP
    _PRIMARY_VARIABLE = "Ta"

    VARIABLE_NAMES = [
        "Ta",
        "RH",
        "water",
        "cloud"
    ]

    def __init__(self, L3G_MET_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L3G_MET_filename)

class L3GSEB(ECOSTRESSGriddedGranule, L3SEB):
    _TILE_CLASS = L3TSEB
    _GRID_NAME = "ECO_L3G_SEB_70m"
    _GRANULE_PREVIEW_CMAP = RN_COLORMAP
    _PRIMARY_VARIABLE = "Rn"

    VARIABLE_NAMES = [
        "Rn",
        "Rg",
        "water",
        "cloud"
    ]

    def __init__(self, L3G_SEB_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L3G_SEB_filename)

class L3GSM(ECOSTRESSGriddedGranule, L3SM):
    _TILE_CLASS = L3TSEB
    _GRID_NAME = "ECO_L3G_SM_70m"
    _GRANULE_PREVIEW_CMAP = SM_COLORMAP
    _PRIMARY_VARIABLE = "SM"

    VARIABLE_NAMES = [
        "SM",
        "water",
        "cloud"
    ]

    def __init__(self, L3G_SM_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L3G_SM_filename)
