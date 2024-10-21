
from abc import ABC, abstractmethod
from datetime import date

from rasters import Raster, RasterGeometry

class VIIRSDownloaderAlbedo(ABC):
    @abstractmethod
    def albedo(
            self,
            date_UTC: date or str,
            geometry: RasterGeometry,
            filename: str = None,
            resampling: str = None) -> Raster:
        pass


class VIIRSDownloaderNDVI(ABC):
    @abstractmethod
    def NDVI(
            self,
            date_UTC: date or str,
            geometry: RasterGeometry,
            filename: str = None,
            resampling: str = None) -> Raster:
        pass
