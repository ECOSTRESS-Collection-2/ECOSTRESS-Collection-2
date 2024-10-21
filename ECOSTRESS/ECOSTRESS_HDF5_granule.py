import logging
import warnings
from abc import abstractmethod
from typing import Any, List

import colored_logging as cl
import h5py
import numpy as np

from ECOSTRESS.ECOSTRESS_granule import ECOSTRESSGranule
from HDFEOS5.HDF5 import HDF5
from HDFEOS5.HDFEOS5 import HDFEOS5
from rasters import Raster, RasterGeometry, KDTree

logger = logging.getLogger(__name__)


class ECOSTRESSHDF5Granule(ECOSTRESSGranule):
    _DATA_GROUP = None
    _PRODUCT_METADATA_GROUP = None
    _STANDARD_METADATA_GROUP = "StandardMetadata"
    _BOUNDARY_WKT_NAME = "SceneBoundaryLatLonWKT"
    _DEFAULT_HDF5_COMPRESSION = "gzip"

    def __init__(
            self,
            product_filename: str,
            data_group_name: str = None,
            *args,
            **kwargs):
        ECOSTRESSGranule.__init__(self, *args, product_filename=product_filename, **kwargs)

        if data_group_name is None:
            data_group_name = self._DATA_GROUP

        self._data_group_name = data_group_name

    @classmethod
    def _clean_hdf_encoding(cls, value):
        if isinstance(value, np.ndarray):
            return [cls._clean_hdf_encoding(item) for item in list(value)]

        if "numpy" in str(type(value)) and hasattr(value, "item"):
            value = value.item()

        if "int" in str(type(value)):
            value = int(value)

        if "float" in str(type(value)):
            value = float(value)

        return value

    @property
    def data_group_name(self):
        if self._data_group_name is None:
            raise ValueError("no data group name given")

        return self._data_group_name

    def read(self, dataset_name: str):
        try:
            with HDF5(self.product_filename, 'r') as file:
                data = file[dataset_name][()]

                if hasattr(data, "decode"):
                    data = str(data.decode())
                elif isinstance(data, np.ndarray) and len(data) == 1:
                    data = data[0]

                return data

        except Exception as e:
            raise IOError(f"unable to read {dataset_name} from: {self.product_filename}")

    @property
    @abstractmethod
    def cloud(self) -> Raster:
        pass

    def data(
            self,
            dataset_name: str,
            apply_scale: bool = True,
            apply_cloud: bool = True) -> Raster:
        if not isinstance(dataset_name, str):
            raise ValueError("invalid dataset name")

        try:
            with HDF5(self.product_filename, 'r') as f:
                data_array = np.array(f[dataset_name])
        except Exception as e:
            raise IOError(f"unable to read dataset {dataset_name} from file {self.product_filename}")

        logger.info(f"dataset {cl.name(dataset_name)} read with dtype {cl.val(data_array.dtype)}")
        variable_name = dataset_name.split('/')[-1]

        if apply_scale:
            dataset_metadata = self.metadata(variable_name)
            
            if "name" in dataset_metadata:
                name = dataset_metadata['name']
            else:
                name = ""
            
            units = dataset_metadata['units']
            fill_value = dataset_metadata['nodata']
            scale = dataset_metadata['scale']
            offset = dataset_metadata['offset']

            logger.debug("name: {}".format(name))
            logger.debug("units: {}".format(units))
            logger.debug("fill: {}".format(fill_value))
            logger.debug("scale: {}".format(scale))
            logger.debug("offset: {}".format(offset))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                raw_min = np.nanmin(data_array)
                raw_max = np.nanmax(data_array)

            logger.debug(f"raw min: {raw_min:0.2f} max: {raw_max:0.2f}")

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                data_array = np.where(data_array == fill_value, np.nan, data_array * scale + offset)
                scaled_min = np.nanmin(data_array)
                scaled_max = np.nanmax(data_array)

            logger.debug(f"scaled min: {scaled_min:0.2f} max: {scaled_max:0.2f}")

            if variable_name.startswith("radiance_"):
                data_array = np.where(data_array < 0.0, np.nan, data_array)

        if apply_cloud and self.cloud is not None:
            data_array = np.where(self.cloud, np.nan, data_array)

        data_raster = Raster(data_array, geometry=self.geometry)

        return data_raster

    def variable(
            self,
            variable_name: str,
            apply_scale: bool = True,
            apply_cloud: bool = True,
            geometry: RasterGeometry = None,
            fill_value: Any = None,
            kd_tree: KDTree = None,
            cmap=None,
            **kwargs):
        if self.data_group_name is None:
            raise ValueError("no data group name")

        data = self.data(
            dataset_name=f"{self.data_group_name}/{variable_name}",
            apply_scale=apply_scale,
            apply_cloud=apply_cloud
        )

        if geometry is not None:
            gridded_data = data.resample(target_geometry=geometry, nodata=fill_value, kd_tree=kd_tree)
            data = gridded_data.astype(data.dtype)

        if cmap is not None:
            print(f"assigning cmap: {cmap}")
            data.cmap = cmap

        return data

    @property
    def variables(self) -> List[str]:
        with HDF5(self.product_filename, "r") as file:
            return file.listing(self.data_group_name)

    def read_metadata_dataset(self, dataset_name):
        return self._clean_hdf_encoding(self.read(dataset_name))

    def metadata(self, variable_name: str) -> dict:
        with HDFEOS5(self.product_filename, "r") as file:
            return file.metadata(f"{self.data_group_name}/{variable_name}")

    @property
    def product_metadata_group_path(self) -> str:
        return self.product_metadata_group_name

    @property
    def standard_metadata_group_path(self) -> str:
        return self.standard_metadata_group_name

    @property
    def standard_metadata(self):
        standard_metadata = {}

        with HDF5(self.product_filename, "r") as file:
            for dataset_name in sorted(file[self.standard_metadata_group_path].keys()):
                full_dataset_name = f"{self.standard_metadata_group_path}/{dataset_name}"
                standard_metadata[dataset_name] = self.read_metadata_dataset(full_dataset_name)

        return standard_metadata

    @property
    def standard_metadata_strings(self):
        return dict([(str(key), str(value)) for key, value in self.standard_metadata.items()])

    @property
    def product_metadata(self):
        product_metadata = {}

        with HDF5(self.product_filename, "r") as file:
            for dataset_name in sorted(file[self.product_metadata_group_path].keys()):
                full_dataset_name = f"{self.product_metadata_group_path}/{dataset_name}"
                product_metadata[dataset_name] = self.read_metadata_dataset(full_dataset_name)

        return product_metadata

    @property
    def product_metadata_strings(self):
        return dict([(str(key), str(value)) for key, value in self.product_metadata.items()])

    def inspect(self):
        with HDF5(self.product_filename, "r") as file:
            file.inspect()

    def copy_metadata(
            self,
            output: h5py.File,
            output_product_metadata_group_name: str = None,
            input_product_metadata_group_name: str = None):
        with h5py.File(self.product_filename, "r") as source:
            dataset_names = []
            source[self._STANDARD_METADATA_GROUP].visit(dataset_names.append)

            standard_metadata = output.create_group(
                f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.standard_metadata_group_name}")

            for dataset_name in dataset_names:
                if "AutomaticQualityFlagExplanation" in dataset_name:
                    continue

                standard_metadata.create_dataset(
                    dataset_name,
                    data=source[f"{self._STANDARD_METADATA_GROUP}/{dataset_name}"]
                )

            if output_product_metadata_group_name is None:
                output_product_metadata_group_name = self.product_metadata_group_name

            if input_product_metadata_group_name is None:
                input_product_metadata_group_name = self.product_metadata_group_name

            dataset_names = []
            source[input_product_metadata_group_name].visit(dataset_names.append)
            product_metadata = output.create_group(f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{output_product_metadata_group_name}")

            for dataset_name in dataset_names:
                if "AutomaticQualityFlagExplanation" in dataset_name:
                    continue

                product_metadata.create_dataset(
                    dataset_name,
                    data=source[f"{input_product_metadata_group_name}/{dataset_name}"]
                )
