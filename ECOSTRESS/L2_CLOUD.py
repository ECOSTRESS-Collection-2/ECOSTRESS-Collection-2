from __future__ import annotations

import logging
import os
from os.path import exists, basename, splitext, abspath, expanduser
from typing import Union, Any, List

import PIL
import h5py
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import Colormap

from ECOSTRESS.ECOSTRESS_granule import ECOSTRESSGranule
from shapely.geometry import Polygon

from ECOSTRESS.ECOSTRESS_HDF5_granule import ECOSTRESSHDF5Granule
from ECOSTRESS.L1B_GEO import L1BGEO
from ECOSTRESS.exit_codes import InputFilesInaccessible
from ECOSTRESS.filenames import find_corresponding_filename
from ECOSTRESS.scan_resampling import resample_scan_by_scan, clip_tails
from rasters import Raster, RasterGeolocation, RasterGeometry, KDTree

# PRIMARY_VARIABLE = "Cloud_test_final"
# PRODUCT_METADATA_GROUP = "L2 CLOUD Metadata"
# GRID_NAME = "ECO_L2G_CLOUD_70m"
#
# L2G_CLOUD_SHORT_NAME = "ECO_L2G_CLOUD"
# L2G_CLOUD_LONG_NAME = "ECOSTRESS Gridded Cloud Mask Instantaneous L2 Global 70 m"
#
# L2T_CLOUD_SHORT_NAME = "ECO_L2T_CLOUD"
# L2T_CLOUD_LONG_NAME = "ECOSTRESS Tiled Cloud Mask Instantaneous L2 Global 70 m"
#
# VARIABLE_NAMES = [
#     "Cloud_test_1",
#     "Cloud_test_2",
#     "Cloud_test_final"
# ]
from timer import Timer

logger = logging.getLogger(__name__)

class L2CLOUDGranule(ECOSTRESSGranule):
    def __init__(self, product_filename: str):
        ECOSTRESSGranule.__init__(self, product_filename=product_filename)

class L2CLOUD(ECOSTRESSHDF5Granule, L2CLOUDGranule):
    _LEVEL = "L2"
    _L2_CLOUD_DATA_GROUP = "SDS"
    _CLOUD_NAME = None
    _LAT_NAME = "latitude"
    _LON_NAME = "longitude"
    _DATA_GROUP = _L2_CLOUD_DATA_GROUP
    _PRODUCT_METADATA_GROUP = "L2 CLOUD Metadata"

    def __init__(
            self,
            L2_CLOUD_filename: str,
            L1B_GEO: L1BGEO = None,
            L1B_GEO_filename: str = None):
        if L2_CLOUD_filename is None:
            raise ValueError("L2 LSTE filename not given")
        elif not exists(L2_CLOUD_filename):
            raise InputFilesInaccessible(f"L2 LSTE file not found: {L2_CLOUD_filename}")

        super(L2CLOUD, self).__init__(product_filename=L2_CLOUD_filename)

        if not isinstance(L2_CLOUD_filename, str):
            raise ValueError("invalid ECOSTRESS L2 CLOUD filename passed when attempting to read geolocation")

        if not exists(L2_CLOUD_filename):
            raise InputFilesInaccessible(f"L2 CLOUD product does not exist: {self.L2_CLOUD_filename}")

        if not os.access(L2_CLOUD_filename, os.R_OK):
            raise InputFilesInaccessible(f"L2 CLOUD product does not have read permission: {self.L2_CLOUD_filename}")

        if L1B_GEO is None:
            if L1B_GEO_filename is None:
                L1B_GEO_filename = find_corresponding_filename(L2_CLOUD_filename, "L1B_GEO")

            L1B_GEO_filename = abspath(expanduser(L1B_GEO_filename))
            L1B_GEO = L1BGEO(L1B_GEO_filename=L1B_GEO_filename)

        self.L1B_GEO = L1B_GEO
        self.L2_CLOUD_filename = L2_CLOUD_filename
        self._cloud = None

    def __repr__(self) -> str:
        return f'L2_CLOUD("{self.L2_CLOUD_filename}")'

    @property
    def Cloud_test_1(self) -> Raster:
        return self.variable("Cloud_test_1", apply_scale=False, apply_cloud=False) == 1

    @property
    def Cloud_test_2(self) -> Raster:
        return self.variable("Cloud_test_2", apply_scale=False, apply_cloud=False) == 1

    @property
    def Cloud_test_final(self) -> Raster:
        return self.variable("Cloud_test_final", apply_scale=False, apply_cloud=False) == 1

    @property
    def Cloud_final(self) -> Raster:
        return self.variable("Cloud_final", apply_scale=False, apply_cloud=False) == 1

    @property
    def Cloud_confidence(self) -> Raster:
        return self.variable("Cloud_confidence", apply_scale=False, apply_cloud=False) == 1

    @property
    def cloud(self) -> Raster:
        if self._cloud is not None:
            return self._cloud

        if "Cloud_final" in self.variables:
            cloud_variable = "Cloud_final"
        elif "Cloud_test_final" in self.variables:
            cloud_variable = "Cloud_test_final"
        else:
            raise IOError("cloud mask dataset not found")

        image = self.variable(cloud_variable, apply_scale=False, apply_cloud=False) == 1
        self._cloud = image

        return image

    def get_browse_image(
            self,
            cmap: Union[Colormap, str] = None,
            mode: str = "RGB") -> PIL.Image.Image:

        if cmap is None:
            cmap = self.granule_preview_cmap

        if isinstance(cmap, str):
            cmap = get_cmap(cmap)

        if "Cloud_final" in self.variables:
            cloud_variable = "Cloud_final"
        elif "Cloud_test_final" in self.variables:
            cloud_variable = "Cloud_test_final"
        else:
            raise IOError("cloud mask dataset not found")

        image = self.variable(cloud_variable)

        browse_image = image.percentilecut.resize(self.granule_preview_shape, resampling="nearest").to_pillow(
            cmap=cmap,
            mode=mode
        )

        return browse_image

    @property
    def boundary(self) -> Polygon:
        return self.geometry.corner_polygon

    @classmethod
    def open(cls, L2_CLOUD_filename: str, L1B_GEO: L1BGEO = None, L1B_GEO_filename: str = None) -> L2CLOUD:
        filename_base = splitext(basename(L2_CLOUD_filename))[0]
        level = filename_base.split("_")[1]
        product = filename_base.split("_")[2]
        build = filename_base.split("_")[-2]

        if level != "L2":
            raise ValueError(f"invalid level in L2 CLOUD filename: {level}")

        if product != "CLOUD":
            raise ValueError(f"invalid product name in L2 CLOUD filename: {product}")

        major_build = int(build[1])

        if major_build == 6:
            return ECOv001L2CLOUD(
                L2_CLOUD_filename=L2_CLOUD_filename,
                L1B_GEO=L1B_GEO,
                L1B_GEO_filename=L1B_GEO_filename
            )
        elif major_build == 7:
            return ECOv002L2CLOUD(
                L2_CLOUD_filename=L2_CLOUD_filename,
                L1B_GEO=L1B_GEO,
                L1B_GEO_filename=L1B_GEO_filename
            )
        else:
            raise ValueError(f"invalid build in L2 CLOUD filename: {build}")

    @property
    def geometry(self) -> RasterGeolocation:
        return self.L1B_GEO.geometry

    def variable(
            self,
            variable_name: str,
            apply_scale: bool = True,
            apply_cloud: bool = True,
            geometry: RasterGeometry = None,
            cell_size_degrees: float = 0.0007,
            search_radius_meters: float = 100,
            nodata: Any = None,
            kd_tree: KDTree = None,
            scan_kd_trees: List[KDTree] = None,
            overlap_strategy: str = "checkerboard",
            **kwargs):
        if self.data_group_name is None:
            raise ValueError("no data group name")

        data = self.data(
            dataset_name=f"{self.data_group_name}/{variable_name}",
            apply_scale=False,
            apply_cloud=False
        )

        dtype = data.dtype

        if nodata is None:
            if "uint8" in str(dtype):
                # nodata = 255
                nodata = 0
            elif "uint16" in str(dtype):
                nodata = 65535
            else:
                raise ValueError(f"unrecognized dtype {dtype}")

        if geometry is not None:
            if overlap_strategy == "checkerboard":
                if kd_tree is None:
                    logger.warning("checkerboard K-D tree not supplied")

                logger.info(f"started resampling {variable_name} with checkerboard strategy")
                timer = Timer()
                gridded_data = data.resample(target_geometry=geometry, nodata=nodata, kd_tree=kd_tree)
                logger.info(f"finished resampling {variable_name} with checkerboard strategy ({timer})")
            elif overlap_strategy == "scan_by_scan":
                if scan_kd_trees is None:
                    logger.warning("scan-by-scan K-D trees not supplied")

                logger.info(f"started resampling {variable_name} with scan-by-scan strategy")
                timer = Timer()

                gridded_data = resample_scan_by_scan(
                    data,
                    target_geometry=geometry,
                    scan_kd_trees=scan_kd_trees,
                    cell_size_degrees=cell_size_degrees,
                    search_radius_meters=search_radius_meters,
                    nodata=nodata
                )

                logger.info(f"finished resampling {variable_name} with scan-by-scan strategy ({timer})")
            elif overlap_strategy == "remove_105_128":
                if kd_tree is None:
                    logger.warning("clipped-tails K-D tree not supplied")

                logger.info(f"started resampling {variable_name} with clipped-tails strategy")
                timer = Timer()
                gridded_data = clip_tails(data).to_geometry(geometry, kd_tree=kd_tree)
                logger.info(f"finished resampling {variable_name} with clipped-tails strategy ({timer})")
            else:
                raise ValueError(f"unrecognized overlap strategy: {overlap_strategy}")

            data = gridded_data.astype(dtype)
            data.nodata = nodata

        return data


class ECOv001L2CLOUD(L2CLOUD):
    _CLOUD_NAME = "CloudMask"

    def __init__(
            self,
            L2_CLOUD_filename: str,
            L1B_GEO: L1BGEO = None,
            L1B_GEO_filename: str = None):
        super(ECOv001L2CLOUD, self).__init__(
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO=L1B_GEO,
            L1B_GEO_filename=L1B_GEO_filename
        )

        self._water = None

    @property
    def CloudMask(self) -> Raster:
        """
        This function reads the QA flag bits of an ECOSTRES L2 CLOUD product as a Raster object.
        :return: <Raster> cloud mask as Raster object
        """
        with h5py.File(self.L2_CLOUD_filename, 'r') as f:
            data = np.array(f[self._L2_CLOUD_DATA_GROUP][self._CLOUD_NAME])

        image = Raster(data, geometry=self.geometry)

        return image

    @property
    def cloud(self) -> Raster:
        """
        This function reads the cloud mask from the QA flag bits of an ECOSTRES L2 CLOUD product as a Raster object.
        :return: <Raster> cloud mask as Raster object
        """
        if self._cloud is not None:
            return self._cloud

        with h5py.File(self.L2_CLOUD_filename, 'r') as f:
            data = (np.array(f[self._L2_CLOUD_DATA_GROUP][self._CLOUD_NAME]).astype(np.uint32) >> 1) & 1 == 1

        image = Raster(data, geometry=self.geometry)
        self._cloud = image

        return image

    @property
    def water(self) -> Raster:
        """
        This function reads the water mask from the QA flag bits of an ECOSTRES L2 CLOUD product as a Raster object.
        :return: <Raster> cloud mask as Raster object
        """
        if self._water is not None:
            return self._water

        if not isinstance(self.L2_CLOUD_filename, str):
            raise ValueError("invalid ECOSTRESS L2 CLOUD filename passed when attempting to read geolocation")

        if not exists(self.L2_CLOUD_filename):
            raise InputFilesInaccessible(f"L2 CLOUD product does not exist: {self.L2_CLOUD_filename}")

        if not os.access(self.L2_CLOUD_filename, os.R_OK):
            raise InputFilesInaccessible(f"L2 CLOUD product does not have read permission: {self.L2_CLOUD_filename}")

        with h5py.File(self.L2_CLOUD_filename, 'r') as f:
            data = (np.array(f[self._L2_CLOUD_DATA_GROUP][self._CLOUD_NAME]).astype(np.uint32) >> 5) & 1 == 1

        image = Raster(data, geometry=self.geometry)

        self._water = image

        return image


class ECOv002L2CLOUD(L2CLOUD):
    # _CLOUD_NAME = "Cloud_test_final"

    @property
    def Cloud_test_1(self) -> Raster:
        """
        This function reads the first cloud test from the ECOSTRESS Collection 2 L2 CLOUD product as a Raster object.
        :return: <Raster> cloud mask as Raster object
        """
        CLOUD_TEST_NAME = "Cloud_test_1"

        if self._cloud is not None:
            return self._cloud

        with h5py.File(self.L2_CLOUD_filename, 'r') as f:
            data = np.array(f[self._L2_CLOUD_DATA_GROUP][CLOUD_TEST_NAME]).astype(np.uint32) == 1

        image = Raster(data, geometry=self.geometry)

        return image

    @property
    def Cloud_test_2(self) -> Raster:
        """
        This function reads the second cloud test from the ECOSTRESS Collection 2 L2 CLOUD product as a Raster object.
        :return: <Raster> cloud mask as Raster object
        """
        CLOUD_TEST_NAME = "Cloud_test_2"

        if self._cloud is not None:
            return self._cloud

        with h5py.File(self.L2_CLOUD_filename, 'r') as f:
            data = np.array(f[self._L2_CLOUD_DATA_GROUP][CLOUD_TEST_NAME]).astype(np.uint32) == 1

        image = Raster(data, geometry=self.geometry)

        return image

    @property
    def Cloud_test_final(self) -> Raster:
        """
        This function reads the cloud mask from the ECOSTRESS Collection 2 L2 CLOUD product as a Raster object.
        :return: <Raster> cloud mask as Raster object
        """
        if self._cloud is not None:
            return self._cloud

        with h5py.File(self.L2_CLOUD_filename, 'r') as f:
            data = np.array(f[self._L2_CLOUD_DATA_GROUP]["Cloud_test_final"]).astype(np.uint32) == 1

        image = Raster(data, geometry=self.geometry)

        return image

    @property
    def Cloud_final(self) -> Raster:
        """
        This function reads the cloud mask from the ECOSTRESS Collection 2 L2 CLOUD product as a Raster object.
        :return: <Raster> cloud mask as Raster object
        """
        if self._cloud is not None:
            return self._cloud

        with h5py.File(self.L2_CLOUD_filename, 'r') as f:
            data = np.array(f[self._L2_CLOUD_DATA_GROUP]["Cloud_final"]).astype(np.uint32) == 1

        image = Raster(data, geometry=self.geometry)

        return image

    @property
    def cloud(self) -> Raster:
        if "Cloud_final" in self.variables:
            return self.Cloud_final
        elif "Cloud_test_final" in self.variables:
            return self.Cloud_test_final
        else:
            raise IOError("cloud mask dataset not found")
