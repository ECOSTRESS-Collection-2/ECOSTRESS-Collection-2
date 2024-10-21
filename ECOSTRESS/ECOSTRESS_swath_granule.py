import logging
from os.path import abspath, expanduser
from typing import List

import h5py
from shapely.geometry import Polygon

from ECOSTRESS.ECOSTRESS_HDF5_granule import ECOSTRESSHDF5Granule
from ECOSTRESS.L1B_GEO import L1BGEO
from ECOSTRESS.L2_CLOUD import L2CLOUD
from ECOSTRESS.filenames import find_corresponding_filename
from rasters import Raster, RasterGeolocation

logger = logging.getLogger(__name__)


class ECOSTRESSSwathGranule(ECOSTRESSHDF5Granule):
    _GRID_NAME = None

    def __init__(
            self,
            product_filename: str,
            L2_CLOUD: L2CLOUD = None,
            L2_CLOUD_filename: str = None,
            L1B_GEO: L1BGEO = None,
            L1B_GEO_filename: str = None,
            data_group_name: str = None,
            **kwargs):
        product_filename = abspath(expanduser(product_filename))

        if L1B_GEO is None:
            if L1B_GEO_filename is None:
                L1B_GEO_filename = find_corresponding_filename(product_filename, "L1B_GEO")

            L1B_GEO_filename = abspath(expanduser(L1B_GEO_filename))
            L1B_GEO = L1BGEO(L1B_GEO_filename=L1B_GEO_filename)

        if L2_CLOUD is None:
            if L2_CLOUD_filename is None:
                L2_CLOUD_filename = find_corresponding_filename(product_filename, "L2_CLOUD")

            L2_CLOUD_filename = abspath(expanduser(L2_CLOUD_filename))
            L2_CLOUD = L2CLOUD.open(L2_CLOUD_filename=L2_CLOUD_filename, L1B_GEO_filename=L1B_GEO_filename)

        ECOSTRESSHDF5Granule.__init__(
            self,
            product_filename=product_filename,
            data_group_name=data_group_name
        )

        self.L1B_GEO = L1B_GEO
        self.L2_CLOUD = L2_CLOUD

    @property
    def L1B_GEO_filename(self) -> str:
        return self.L1B_GEO.L1B_GEO_filename

    # @property
    # def L2_LSTE_filename(self) -> str:
    #     return self.product_filename

    @property
    def L2_CLOUD_filename(self) -> str:
        return self.L2_CLOUD.L2_CLOUD_filename

    @property
    def cloud(self) -> Raster:
        return self.L2_CLOUD.cloud

    @property
    def water(self) -> Raster:
        return self.L2_CLOUD.water

    @property
    def filename(self) -> str:
        return self.product_filename

    @property
    def gridded_level(self) -> str:
        return f"{self.level}G"

    @property
    def gridded_granule_name(self) -> str:
        granule_name = self.granule_name
        granule_name = granule_name.replace("ECOSTRESS_", "ECOv002_")
        granule_name = granule_name.replace(f"_{self.level}_", f"_{self.gridded_level}_")

        return granule_name

    @property
    def cloud_filename(self) -> str:
        return self.L2_CLOUD.L2_CLOUD_filename

    def visit(self, dataset: str = None) -> List[str]:
        result = []

        with h5py.File(self.product_filename, "r") as file:
            if dataset is None:
                file.visit(result.append)
            else:
                file[dataset].visit(result.append)

        return result

    @property
    def geometry(self) -> RasterGeolocation:
        return self.L1B_GEO.geometry

    @property
    def boundary(self) -> Polygon:
        return self.geometry.corner_polygon
