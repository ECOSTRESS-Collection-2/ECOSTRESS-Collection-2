from datetime import datetime
from typing import Union

from dateutil import parser

from ECOSTRESS_colors import NDVI_COLORMAP
from .ECOSTRESS_granule import ECOSTRESSGranule
from .ECOSTRESS_gridded_tiled_granule import ECOSTRESSTiledGranule

__author__ = "Gregory Halverson"

PRIMARY_VARIABLE = "NDVI"
PREVIEW_CMAP = NDVI_COLORMAP


class InaccessibleL2STARSError(Exception):
    pass


class L2STARSGranule(ECOSTRESSGranule):
    _PRODUCT_NAME = "L2T_STARS"
    _PRIMARY_VARIABLE = "NDVI"
    _GRANULE_PREVIEW_CMAP = PREVIEW_CMAP

    def __init__(self, product_filename: str):
        super(L2STARSGranule, self).__init__(product_filename=product_filename)

        self._NDVI = None
        self._NDVI_UQ = None
        self._albedo = None
        self._albedo_UQ = None

    @property
    def NDVI(self):
        if self._NDVI is None:
            self._NDVI = self.variable("NDVI")

        return self._NDVI

    @property
    def NDVI_UQ(self):
        if self._NDVI_UQ is None:
            self._NDVI_UQ = self.variable("NDVI-UQ")

        return self._NDVI_UQ

    @property
    def albedo(self):
        if self._albedo is None:
            self._albedo = self.variable("albedo")

        return self._albedo

    @property
    def albedo_UQ(self):
        if self._albedo_UQ is None:
            self._albedo_UQ = self.variable("albedo-UQ")

        return self._albedo_UQ

    @property
    def orbit(self):
        return None

    @property
    def scene(self):
        return None


class L2TSTARS(ECOSTRESSTiledGranule, L2STARSGranule):
    _PRODUCT_NAME = "L2T_STARS"
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE
    _GRANULE_PREVIEW_CMAP = PREVIEW_CMAP

    def __init__(
            self,
            product_location: str = None,
            orbit: int = None,
            scene: int = None,
            tile: str = None,
            time_UTC: Union[datetime, str] = None,
            build: str = None,
            process_count: int = None,
            containing_directory: str = None):
        L2STARSGranule.__init__(self, product_filename=product_location)

        ECOSTRESSTiledGranule.__init__(
            self,
            orbit=orbit,
            scene=scene,
            tile=tile,
            time_UTC=time_UTC,
            build=build,
            process_count=process_count,
            product_location=product_location,
            containing_directory=containing_directory
        )

    @property
    def orbit(self):
        return None

    @property
    def scene(self):
        return None

    @classmethod
    def generate_granule_name(
            cls,
            product_name: str,
            orbit: int,
            scene: int,
            tile: str,
            time_UTC: Union[datetime, str],
            build: str,
            process_count: int):
        if product_name is None:
            raise ValueError("invalid product name")

        if orbit is None:
            raise ValueError("invalid orbit")

        if scene is None:
            raise ValueError("invalid scene")

        if tile is None:
            raise ValueError("invalid tile")

        if time_UTC is None:
            raise ValueError("invalid time")

        if build is None:
            raise ValueError("invalid build")

        if process_count is None:
            raise ValueError("invalid process count")

        if isinstance(time_UTC, str):
            time_UTC = parser.parse(time_UTC)

        granule_name = f"ECOv002_{product_name}_{tile}_{time_UTC:%Y%m%dT%H%M%S}_{build}_{process_count:02d}"

        return granule_name
