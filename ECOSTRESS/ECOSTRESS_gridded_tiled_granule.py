from __future__ import annotations

import json
import logging
import os
import posixpath
import shutil
import zipfile
from abc import abstractmethod
from datetime import datetime, timedelta
from glob import glob
from os import makedirs
from os.path import dirname, join, exists, splitext, basename, abspath, expanduser
from time import process_time, perf_counter
from typing import List, Union, Optional

import geopandas as gpd
import h5py
import numpy as np
from ECOSTRESS.exit_codes import BlankOutput
from dateutil import parser
from matplotlib.cm import get_cmap
from matplotlib.colors import Colormap
from shapely import wkt
from shapely.geometry import Polygon

import colored_logging as cl
import rasters
import rasters as rt
from ECOSTRESS.ECOSTRESS_HDF5_granule import ECOSTRESSHDF5Granule
from ECOSTRESS.ECOSTRESS_granule import ECOSTRESSGranule, DEFAULT_JSON_INDENT
from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule
from ECOSTRESS.XML_metadata import write_XML_metadata

from HDFEOS5.HDFEOS5 import HDFEOS5, h5py_copy
from he5py import he5py
from rasters import RasterGrid, Raster, RasterGeometry
from sentinel_tile_grid import SentinelTileGrid, sentinel_tile_grid
from timer import Timer

logger = logging.getLogger(__name__)

class ECOSTRESSTiledGranule(ECOSTRESSGranule):
    _PRODUCT_NAME = None

    _COMPRESSION = "zstd"
    _GRANULE_PREVIEW_CMAP = "jet"

    _STANDARD_METADATA_GROUP = "StandardMetadata"
    _PRODUCT_METADATA_GROUP = "ProductMetadata"

    def __init__(
            self,
            product_location: str = None,
            orbit: int = None,
            scene: int = None,
            tile: str = None,
            time_UTC: Union[datetime, str] = None,
            build: str = None,
            process_count: int = None,
            *args,
            product_name: str = None,
            containing_directory: str = None,
            compression: str = None,
            layer_preview_quality: int = None,
            granule_preview_cmap: Union[Colormap, str] = None,
            granule_preview_shape: (int, int) = None,
            granule_preview_quality: int = None,
            **kwargs):
        ECOSTRESSGranule.__init__(
            self,
            *args,
            product_filename=product_location,
            granule_preview_cmap=granule_preview_cmap,
            granule_preview_shape=granule_preview_shape,
            granule_preview_quality=granule_preview_quality,
            **kwargs
        )

        if product_name is None:
            product_name = self._PRODUCT_NAME

        if product_name is None:
            raise ValueError("product name not given for ECOSTRESS tiled product")

        if product_location is None:
            if isinstance(time_UTC, str):
                time_UTC = parser.parse(time_UTC)

            granule_name = self.generate_granule_name(
                product_name=product_name,
                orbit=orbit,
                scene=scene,
                tile=tile,
                time_UTC=time_UTC,
                build=build,
                process_count=process_count
            )

            if containing_directory is None:
                containing_directory = "."

            product_location = join(containing_directory, granule_name)
        else:
            granule_name = splitext(basename(product_location))[0]

        product_location = abspath(expanduser(product_location))

        if not exists(product_location) and not product_location.endswith(".zip"):
            makedirs(product_location, exist_ok=True)

        if layer_preview_quality is None:
            layer_preview_quality = self._LAYER_PREVIEW_QUALITY

        if compression is None:
            compression = self._COMPRESSION

        self._metadata_dict = None

        self._product_location = product_location
        self._orbit = orbit
        self._scene = scene
        self._tile = tile
        self._time_UTC = time_UTC
        self._time_solar = None
        self._build = build
        self._process_count = process_count
        self._granule_name = granule_name
        self.layer_preview_quality = layer_preview_quality
        self.compression = compression

    @property
    def product_location(self) -> str:
        return self._product_location

    @property
    def filename(self) -> str:
        return self.product_location

    @property
    def product_directory(self) -> str:
        product_location = self.product_location

        if self.is_zip:
            product_directory = splitext(product_location)[0]
        else:
            product_directory = product_location

        return product_directory

    def __repr__(self):
        return f'{self.__class__.__name__}("{self.product_location}")'

    @property
    def metadata_filename_base(self) -> str:
        return f"{self.granule_name}.json"

    @property
    def metadata_filename(self) -> str:
        return join(self.product_directory, self.metadata_filename_base)

    @property
    def metadata_dict(self) -> dict:
        if self._metadata_dict is not None:
            return self._metadata_dict

        if self.is_zip:
            zip_filename = self.product_location
            filename = posixpath.join(self.granule_name, self.metadata_filename_base)

            logger.info(f"reading metadata file {cl.file(filename)} from product zip {cl.file(zip_filename)}")

            with zipfile.ZipFile(zip_filename) as file:
                JSON_text = file.read(filename)
        else:
            filename = self.metadata_filename

            with open(filename, "r") as file:
                JSON_text = file.read()

        metadata_dict = json.loads(JSON_text)
        self._metadata_dict = metadata_dict

        return metadata_dict

    @property
    def standard_metadata(self) -> dict:
        return self.metadata_dict[self.standard_metadata_group_name]

    @property
    def product_metadata(self) -> dict:
        return self.metadata_dict[self.product_metadata_group_name]

    @property
    def water(self) -> Raster:
        return self.variable("water")

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

        granule_name = f"ECOv002_{product_name}_{orbit:05d}_{scene:03d}_{tile}_{time_UTC:%Y%m%dT%H%M%S}_{build}_{process_count:02d}"

        return granule_name

    def add_layer(
            self,
            variable_name: str,
            image: Raster,
            cmap: Union[Colormap, str] = None,
            include_COG: bool = True,
            include_geojpeg: bool = True) -> str:
        if cmap is None:
            cmap = "jet"

        if isinstance(cmap, str):
            cmap = get_cmap(cmap)

        # print("add_layer")
        filename_base = self.layer_filename(variable_name)
        filename = join(self.product_location, filename_base)
        # print(f"filename: {filename}")
        logger.info(f"adding {cl.name(self.product)} layer {cl.name(variable_name)} min: {cl.val(np.nanmin(image))} mean: {cl.val(np.nanmean(image))} max: {cl.val(np.nanmax(image))} cmap: {cl.name(cmap.name)} file: {cl.file(filename)}")

        preview_filename = filename.replace(".tif", ".jpeg")

        if include_COG:
            image.to_COG(
                filename=filename,
                compress=self.compression,
                preview_filename=preview_filename,
                preview_quality=self.layer_preview_quality,
                cmap=cmap
            )

        if include_COG:
            return filename

        if include_geojpeg:
            return preview_filename

    def write_metadata(self, metadata_dict: dict, indent=DEFAULT_JSON_INDENT):
        filename = self.metadata_filename
        JSON_text = json.dumps(metadata_dict, indent=indent)

        with open(filename, "w") as file:
            file.write(JSON_text)

    @property
    def storage(self):
        return os.stat(self.product_filename).st_size

    def write_zip(self, zip_filename: str):
        product_directory = self.product_directory
        logger.info(f"writing product zip file: {cl.file(product_directory)} -> {cl.file(zip_filename)}")

        directory_name = splitext(basename(zip_filename))[0]

        with zipfile.ZipFile(zip_filename, "w") as zip_file:
            for filename in glob(join(product_directory, "*")):
                arcname = join(directory_name, basename(filename))
                zip_file.write(
                    filename=filename,
                    arcname=arcname
                )

        XML_metadata_filename = zip_filename.replace(".zip", ".zip.xml")
        logger.info(f"writing XML metadata file: {cl.file(XML_metadata_filename)}")
        write_XML_metadata(self.standard_metadata, XML_metadata_filename)

        # TODO should be zip-file validator

        if not exists(zip_filename):
            raise IOError(f"unable to create tiled product zip: {zip_filename}")

    @classmethod
    def from_scene(
            cls,
            gridded_granule: ECOSTRESSGriddedGranule,
            tile: str,
            tile_granule_directory: str = None,
            tile_granule_name: str = None,
            geometry: RasterGeometry = None,
            variables: List[str] = None,
            compression: str = None,
            overwrite: bool = False,
            skip_blank: bool = True) -> ECOSTRESSTiledGranule or None:
        if compression is None:
            compression = cls._COMPRESSION

        if tile_granule_name is None:
            tile_granule_name = gridded_granule.tile_granule_name(tile)

        if tile_granule_directory is None:
            tile_granule_directory = tile_granule_name

        logger.info(
            f"target granule directory: {cl.dir(tile_granule_directory)}"
        )

        metadata = gridded_granule.metadata_dict
        metadata["StandardMetadata"]["LocalGranuleID"] = f"{tile_granule_name}.zip"

        orbit = gridded_granule.orbit
        scene = gridded_granule.scene
        time_UTC = gridded_granule.time_UTC
        build = gridded_granule.build
        process_count = gridded_granule.process_count

        granule = cls(
            product_location=tile_granule_directory,
            orbit=orbit,
            scene=scene,
            tile=tile,
            time_UTC=time_UTC,
            build=build,
            process_count=process_count,
            compression=compression
        )

        granule.write_metadata(metadata)

        if variables is None:
            output_variables = gridded_granule.tiled_output_variables
        else:
            output_variables = variables

        for j, variable in enumerate(output_variables):
            logger.info(f"processing variable: {variable}")
            output_filename = join(tile_granule_directory, f"{tile_granule_name}_{variable}.tif")

            if exists(output_filename) and not overwrite:
                logger.warning(f"file already exists: {cl.file(output_filename)}")
                continue

            logger.info(
                f"started processing variable {variable} ({j + 1} / {len(output_variables)}) "
                f"for granule: {tile_granule_name}"
            )

            timer = Timer()

            image = gridded_granule.variable(
                variable,
                apply_scale=True,
                apply_cloud=True,
                geometry=geometry
            )

            if skip_blank and np.all(np.isnan(image)):
                raise BlankOutput(f"blank output for layer {variable} at tile {tile} at time {time_UTC}")

            granule.add_layer(variable, image)

            logger.info(
                f"finished processing variable {variable} ({j + 1} / {len(output_variables)}) "
                f"for granule: {tile_granule_name} "
                f"({cl.time(timer)})"
            )

        return granule

    @classmethod
    def scan_directory(cls, directory: str) -> List[ECOSTRESSTiledGranule]:
        filenames = glob(join(directory, "*"))
        filenames = filter(lambda filename: splitext(basename(filename))[0].split("_")[1].endswith("T"), filenames)

        granules = []

        for filename in filenames:
            try:
                granule = cls(filename)
                granules.append(granule)
            except:
                continue

        return granules

    @property
    def geometry(self) -> RasterGrid:
        if self._geometry is not None:
            return self._geometry

        URI = self.layer_URI(self.primary_variable)
        grid = RasterGrid.open(URI)

        return grid

    @property
    def boundary(self) -> Polygon:
        return self.geometry.boundary

    @property
    def orbit(self) -> int:
        if self._orbit is not None:
            return self._orbit
        else:
            return int(self.granule_name.split('_')[-6])

    @property
    def scene(self) -> int:
        if self._scene is not None:
            return self._scene
        else:
            return int(self.granule_name.split('_')[-5])

    @property
    def tile(self) -> str:
        if self._tile is not None:
            return self._tile
        else:
            return self.granule_name.split('_')[-4]

    @property
    def time_UTC(self) -> datetime:
        if self._time_UTC is not None:
            return self._time_UTC
        else:
            return parser.parse(self.granule_name.split('_')[-3])

    def UTC_to_solar(self, time_UTC: datetime, lon: float) -> datetime:
        return time_UTC + timedelta(hours=(np.radians(lon) / np.pi * 12))

    @property
    def time_solar(self) -> datetime:
        if self._time_solar is not None:
            return self._time_solar
        else:
            return self.UTC_to_solar(self.time_UTC, self.geometry.centroid_latlon.x)

    @property
    def build(self) -> str:
        if self._build is not None:
            return self._build
        else:
            return self.granule_name.split('_')[-2]

    @property
    def process_count(self) -> int:
        if self._process_count is not None:
            return self._process_count
        else:
            return int(self.granule_name.split('_')[-1])

    @property
    def product(self) -> str:
        return "_".join(self.granule_name.split("_")[1:-6])

    @property
    def granule_name(self) -> str:
        if self._granule_name is not None:
            return self._granule_name
        else:
            return self.generate_granule_name(
                product_name=self.product,
                orbit=self.orbit,
                scene=self.scene,
                tile=self.tile,
                time_UTC=self.time_UTC,
                build=self.build,
                process_count=self.process_count
            )

    @property
    def filenames(self) -> List[str]:
        # return glob(join(self.product_directory, "*.tif"))
        with zipfile.ZipFile(self.product_filename) as zip_file:
            return zip_file.namelist()

    @property
    def layer_filenames(self) -> List[str]:
        return [
            filename
            for filename
            in self.filenames
            if filename.endswith(".tif")
        ]

    def URI_for_filename(self, filename: str) -> str:
        return f"zip://{abspath(expanduser(self.product_filename))}!/{self.granule_name}/{filename}"

    @property
    def layer_URIs(self) -> List[str]:
        return [
            self.URI_for_filename(filename)
            for filename
            in self.layer_filenames
        ]

    @property
    def variables(self) -> List[str]:
        return [
            str(splitext(basename(filename))[0].split("_")[-1])
            for filename
            in self.filenames
            if filename.endswith(".tif")
        ]

    def layer_filename(self, variable: str) -> Optional[str]:
        if not isinstance(variable, str):
            raise ValueError("invalid variable")

        return f"{self.granule_name}_{variable}.tif"

    @property
    def is_zip(self) -> bool:
        return self.product_location.endswith(".zip")

    def layer_URI(self, variable: str) -> Optional[str]:
        # return join(self.product_directory, f"{self.granule_name}_{variable}.tif")

        layer_filename = self.layer_filename(variable)

        if layer_filename is None:
            return None

        if self.is_zip:
            URI = self.URI_for_filename(layer_filename)
        else:
            URI = join(self.product_location, layer_filename)

        return URI

    def variable(self, variable: str, geometry: RasterGeometry = None, cmap=None, **kwargs) -> Raster:
        URI = self.layer_URI(variable)
        logger.info(f"started reading {self.product} {variable}: {cl.URL(URI)}")
        start_time = perf_counter()
        image = Raster.open(URI)

        if geometry is not None:
            logger.info(f"projecting {self.product} {variable}")
            image = image.to_geometry(geometry, **kwargs)

        end_time = perf_counter()
        duration = end_time - start_time
        logger.info(f"finished reading {self.product} {variable} ({duration:0.2f}s)")

        if cmap is not None:
            print(f"assigning cmap: {cmap}")
            image.cmap = cmap

        if "float" in str(image.dtype):
            image.nodata = np.nan

        return image

    @classmethod
    def mosaic(cls, variable: str, granules: List[ECOSTRESSTiledGranule], geometry: RasterGeometry):
        return rasters.mosaic((granule.variable(variable) for granule in granules), geometry)


class ECOSTRESSGriddedGranule(ECOSTRESSHDF5Granule):
    _TILE_CLASS = ECOSTRESSTiledGranule
    _GRID_NAME = None
    _PRIMARY_VARIABLE = None
    _DEFAULT_GEOTIFF_COMPRESSION = "zstd"
    _DEFAULT_COMPRESSION_LEVEL = 17
    _DEFAULT_HDF5_SCALEOFFSET = 2
    _DEFAULT_GRIDDED_PRODUCT_APPLY_SCALE = True
    _DEFAULT_GRIDDED_PRODUCT_APPLY_CLOUD = False
    _DEFAULT_PROJECTION_SYSTEM = "global_geographic"
    _DEFAULT_OVERWRITE = False
    _DATASET_NAME_TRANSLATIONS = {}
    VARIABLE_NAMES = []

    UINT8_FILL = 255
    UINT16_FILL = 65535

    def __init__(
            self,
            product_filename: str,
            target_resolution: float = None):
        ECOSTRESSHDF5Granule.__init__(self, product_filename=product_filename)

        if target_resolution is None:
            target_resolution = self._DEFAULT_UTM_CELL_SIZE

        self.target_resolution = target_resolution
        self.sentinel = SentinelTileGrid(target_resolution=target_resolution)

    def __repr__(self):
        return f'{self.__class__.__name__}("{self.product_filename}")'

    @property
    def tiled_level(self):
        return self.level.replace("G", "T")

    @property
    def tiled_scene_name(self):
        return self.granule_name.replace(self.level, self.tiled_level)

    @property
    def grid_name(self) -> str:
        if self._GRID_NAME is not None:
            return self._GRID_NAME
        else:
            raise ValueError(f"grid name not defined for: {self.__class__.__name__}")

    @property
    def data_group_name(self) -> str:
        return f"HDFEOS/GRIDS/{self.grid_name}/Data Fields"

    @property
    def grid(self) -> RasterGrid:
        with HDFEOS5(self.product_filename, "r") as file:
            return file.grid(self.grid_name)

    @property
    def geometry(self) -> RasterGrid:
        if self._geometry is None:
            self._geometry = self.grid

        return self._geometry

    @property
    def product_metadata_group_path(self) -> str:
        return f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.product_metadata_group_name}"

    @property
    def standard_metadata_group_path(self) -> str:
        return f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.standard_metadata_group_name}"

    @property
    def boundary_WKT_name(self) -> str:
        return self._BOUNDARY_WKT_NAME

    @property
    def boundary_WKT(self) -> str:
        with HDFEOS5(self.product_filename, "r") as file:
            value = file[self.standard_metadata_group_path][self.boundary_WKT_name][()]

            if isinstance(value, str):
                return value
            elif isinstance(value, bytes):
                return value.decode()
            else:
                raise ValueError(f"invalid {self.boundary_WKT_name} type: {type(value)}")

    @property
    def boundary(self) -> Polygon:
        return wkt.loads(self.boundary_WKT)

    @property
    @abstractmethod
    def group_name(self):
        pass

    def metadata(self, variable_name: str):
        with HDFEOS5(self.product_filename, "r") as file:
            return file.metadata(f"{self.data_group_name}/{variable_name}")

    def scale(self, dataset_name: str) -> (float or int, float or int, float or int):
        with HDFEOS5(self.product_filename, "r") as file:
            return file.scale(dataset_name, self.grid_name)

    @property
    @abstractmethod
    def cloud(self) -> Raster:
        pass

    def variable(
            self,
            variable_name: str,
            apply_scale: bool = True,
            apply_cloud: bool = True,
            geometry: RasterGeometry = None,
            **kwargs) -> Raster:

        logger.info(
            f"started reading HDF-EOS5 dataset: {cl.name(variable_name)} "
            f"file: {cl.file(self.product_filename)}"
        )

        HDF5_read_timer = Timer()

        with HDFEOS5(self.product_filename, "r") as file:
            image = file.variable(
                dataset_name=variable_name,
                grid_name=self.grid_name,
                apply_scale=apply_scale,
                geometry=geometry
            )

        logger.info(
            f"finished reading HDFEOS5 dataset: {cl.name(variable_name)} "
            f"file: {cl.file(self.product_filename)} ({cl.time(HDF5_read_timer)})"
        )

        if apply_cloud and "cloud" in self.variables:
            cloud = self.variable(
                variable_name="cloud",
                apply_scale=False,
                apply_cloud=False,
                geometry=geometry
            )

            image = rasters.where(cloud, np.nan, image)

        return image

    @classmethod
    def gridded_filename_base(cls, swath_granule: ECOSTRESSSwathGranule) -> str:
        return f"{swath_granule.gridded_granule_name}.h5"

    @property
    def product_metadata_group_name(self):
        return "ProductMetadata"

    def to_geotiff(self, directory: str = None):
        if directory is None:
            directory = splitext(self.product_filename)[0]

        makedirs(directory, exist_ok=True)

        for variable in self.variables:
            filename = join(directory, f"{self.granule_name}_{variable}.tif")
            image = self.variable(variable)
            logger.info(f"writing gridded layer {cl.name(variable)} to GeoTIFF: {cl.file(filename)}")
            image.to_geotiff(filename)

        metadata = self.metadata_dict
        metadata_filename = join(directory, f"{self.granule_name}.json")

        with open(metadata_filename, "w") as file:
            file.write(json.dumps(metadata, indent=2))

    def copy_metadata(
            self,
            output: h5py.File,
            standard_metadata_group_name: str = None,
            standard_metadata_additional: dict = None,
            standard_metadata_exclude: List[str] = None,
            product_metadata_group_name: str = None,
            product_metadata_additional: dict = None,
            product_metadata_exclude: List[str] = None):
        if standard_metadata_group_name is None:
            standard_metadata_group_name = self.standard_metadata_group_name

        if product_metadata_group_name is None:
            product_metadata_group_name = self.product_metadata_group_name

        with h5py.File(self.product_filename, "r") as source:
            dataset_names = []
            source[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.standard_metadata_group_name}"].visit(dataset_names.append)

            standard_metadata = output.create_group(
                f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.standard_metadata_group_name}")

            if standard_metadata_additional is not None:
                for key in standard_metadata_additional.keys():
                    dataset_names.append(key)

            dataset_names = sorted(set(dataset_names))

            for dataset_name in dataset_names:
                dataset_name = str(dataset_name)

                if dataset_name.startswith("Cloud"):
                    continue

                if dataset_name.startswith("Emis"):
                    continue

                if dataset_name.startswith("LST"):
                    continue

                if dataset_name in standard_metadata_additional.keys():
                    if standard_metadata_exclude is not None and dataset_name in standard_metadata_exclude:
                        continue

                    value = standard_metadata_additional[dataset_name]
                    logger.info(f"adding standard metadata for {dataset_name}: {value}")
                    standard_metadata.create_dataset(
                        dataset_name,
                        data=standard_metadata_additional[dataset_name]
                    )
                else:
                    value = source[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{standard_metadata_group_name}/{dataset_name}"]
                    logger.info(f"adding standard metadata for {dataset_name}: {value}")
                    standard_metadata.create_dataset(
                        dataset_name,
                        data=source[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{standard_metadata_group_name}/{dataset_name}"]
                    )

            # for dataset_name in dataset_names:
            #     standard_metadata.create_dataset(
            #         dataset_name,
            #         data=source[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.standard_metadata_group_name}/{dataset_name}"]
            #     )

            dataset_names = []
            source[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{product_metadata_group_name}"].visit(dataset_names.append)
            product_metadata = output.create_group(f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{product_metadata_group_name}")

            if product_metadata_additional is not None:
                for key in product_metadata_additional.keys():
                    dataset_names.append(key)

            dataset_names = sorted(set(dataset_names))

            for dataset_name in dataset_names:
                dataset_name = str(dataset_name)

                if dataset_name.startswith("Cloud"):
                    continue

                if dataset_name.startswith("Emis"):
                    continue

                if dataset_name.startswith("LST"):
                    continue

                if dataset_name in product_metadata_additional.keys():
                    if product_metadata_exclude is not None and dataset_name in product_metadata_exclude:
                        continue

                    value = product_metadata_additional[dataset_name]
                    logger.info(f"adding product metadata for {dataset_name}: {value}")
                    product_metadata.create_dataset(
                        dataset_name,
                        data=product_metadata_additional[dataset_name]
                    )
                else:
                    value = source[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{product_metadata_group_name}/{dataset_name}"]
                    logger.info(f"adding product metadata for {dataset_name}: {value}")
                    product_metadata.create_dataset(
                        dataset_name,
                        data=source[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{product_metadata_group_name}/{dataset_name}"]
                    )


    @classmethod
    def from_tiles(
            cls,
            filename: str,
            tile_filenames: List[str],
            gridded_source_granule: ECOSTRESSGriddedGranule,
            grid_name: str = None,
            standard_metadata_additional: dict = None,
            product_metadata_additional: dict = None,
            product_metadata_exclude: List[str] = None,
            ancillary_filenames_name: str = None,
            geotiff_diagnostics: bool = False) -> ECOSTRESSGriddedGranule:
        if grid_name is None:
            grid_name = cls._GRID_NAME

        if product_metadata_additional is None:
            product_metadata_additional = {}

        if ancillary_filenames_name is None:
            ancillary_filenames_name = "AncillaryNWP"

        if standard_metadata_additional is None:
            standard_metadata_additional = {}

        geometry = gridded_source_granule.geometry

        logger.info(f"mosaicking {len(tile_filenames)} tiles into composite: {filename}")

        output_filename_partial = f"{filename}.partial"
        output_filename_rewritten = f"{filename}.rewritten"

        if exists(output_filename_partial):
            logger.warning(f"removing prior partial file: {output_filename_partial}")
            os.remove(output_filename_partial)

        makedirs(dirname(output_filename_partial), exist_ok=True)

        ancillary_filenames = set([])

        logger.info(f"opening partial HDF-EOS5 file for writing: {cl.file(output_filename_partial)}")

        with he5py.File(output_filename_partial, "w") as output_file:
            output_grid = output_file.create_grid(grid_name=grid_name, geometry=geometry)

            first_granule = cls._TILE_CLASS(tile_filenames[0])

            for variable in cls.VARIABLE_NAMES:
                logger.info(f"processing {cl.name(variable)} mosaic {geometry.shape} for: {filename}")

                dtype = first_granule.variable(variable).dtype
                logger.info(f"using dtype: {cl.val(dtype)}")

                if "float" in str(dtype):
                    composite_sum = Raster(np.full(geometry.shape, 0, dtype=dtype), geometry=geometry)
                    composite_count = Raster(np.full(geometry.shape, 0, dtype=np.uint16), geometry=geometry)

                    for tile_filename in sorted(tile_filenames):
                        logger.info(f"loading {cl.name(variable)} from tile: {cl.file(tile_filename)}")
                        tile_granule = cls._TILE_CLASS(tile_filename)
                        ancillary_filenames = set(ancillary_filenames) | set(tile_granule.product_metadata[ancillary_filenames_name].split(","))
                        tile_image = tile_granule.variable(variable)
                        logger.info(f"loaded tile image {cl.val(tile_image.shape)}: {cl.file(tile_filename)}")
                        projected_tile_image = tile_image.to_geometry(geometry)
                        composite_sum = rt.where(np.isnan(projected_tile_image), composite_sum,
                                                 composite_sum + projected_tile_image)
                        projected_tile_count = rt.where(np.isnan(projected_tile_image), 0, 1)
                        composite_count = composite_count + projected_tile_count
                        logger.info(
                            f"composite fill: {(np.count_nonzero(composite_count) / composite_count.size * 100):0.2f}%")

                    composite_image = rt.where(composite_count > 0, composite_sum / composite_count, np.nan)

                else:
                    composite_image = Raster(np.full(geometry.shape, 0), geometry=geometry)

                    for tile_filename in sorted(tile_filenames):
                        logger.info(f"loading {cl.name(variable)} from tile: {cl.file(tile_filename)}")
                        tile_granule = cls._TILE_CLASS(tile_filename)
                        tile_image = tile_granule.variable(variable)
                        logger.info(f"loaded tile image {cl.val(tile_image.shape)}: {cl.file(tile_filename)}")
                        projected_tile_image = tile_image.to_geometry(geometry)
                        composite_image = rt.where(np.isnan(projected_tile_image), composite_image,
                                                   projected_tile_image)
                        logger.info(
                            f"composite fill: {(np.count_nonzero(~np.isnan(composite_image)) / composite_image.size * 100):0.2f}%")

                logger.info(
                    f"composite missing: {(np.count_nonzero(np.isnan(composite_image)) / composite_image.size * 100):0.2f}%")
                logger.info(f"writing {cl.name(variable)} composite: {cl.file(output_filename_partial)}")

                output_grid.write_float(
                    field_name=variable,
                    image=composite_image
                )

                if variable == cls._PRIMARY_VARIABLE:
                    percent_good_quality = 100 * (1 - np.count_nonzero(np.isnan(composite_image)) / composite_image.size)
                    product_metadata_additional["QAPercentGoodQuality"] = percent_good_quality
                elif variable == "cloud":
                    percent_cloud = 100 * np.count_nonzero(composite_image) / composite_image.size
                    product_metadata_additional["QAPercentCloudCover"] = percent_cloud

        product_metadata_additional = {
            ancillary_filenames_name: ",".join(ancillary_filenames)
        }

        standard_metadata_additional["LocalGranuleID"] = basename(filename)

        with HDFEOS5(output_filename_partial, "r+") as output_file:
            logger.info(f"started copying ECOSTRESS metadata to HDF-EOS5 file: {cl.file(output_filename_partial)}")
            gridded_source_granule.copy_metadata(
                output=output_file,
                standard_metadata_additional=standard_metadata_additional,
                product_metadata_additional=product_metadata_additional,
                product_metadata_exclude=product_metadata_exclude
            )

        # shutil.move(output_filename_partial, filename)
        logger.info(f"re-writing HDF-EOS5 file: {output_filename_partial} -> {output_filename_rewritten}")
        h5py_copy(output_filename_partial, output_filename_rewritten)
        logger.info(f"removing original HDF-EOS5 output: {output_filename_partial}")
        os.remove(output_filename_partial)
        logger.info(f"renaming completed output file: {output_filename_rewritten} -> {filename}")
        shutil.move(output_filename_rewritten, filename)

        logger.info(f"completed HDF-EOS5 generation: {cl.file(filename)}")
        granule = cls(filename)

        confirmed_output_variables = granule.variables

        logger.info(f"checking output layers in product file: {cl.file(filename)}")

        for variable in cls.VARIABLE_NAMES:
            if variable in confirmed_output_variables:
                logger.info(f"* {cl.name(variable)} confirmed")
            else:
                raise IOError(f"output layer {variable} not found: {filename}")

        if geotiff_diagnostics:
            granule.to_geotiff()

        XML_metadata_filename = filename.replace(".h5", ".h5.met")
        logger.info(f"writing XML metadata file: {cl.file(XML_metadata_filename)}")
        write_XML_metadata(granule.standard_metadata, XML_metadata_filename)

        return granule

    @classmethod
    def translate_dataset_name(cls, dataset_name: str) -> str:
        if dataset_name in cls._DATASET_NAME_TRANSLATIONS:
            return cls._DATASET_NAME_TRANSLATIONS[dataset_name]
        else:
            return dataset_name

    @classmethod
    def from_swath(
            cls,
            swath_granule: ECOSTRESSSwathGranule,
            output_filename: str = None,
            output_directory: str = None,
            projection_system: str = None,
            cell_size: float = None,
            grid_name: str = None,
            compression: str = None,
            scaleoffset: int = None,
            apply_scale: bool = None,
            apply_cloud: bool = None,
            overwrite: bool = None) -> ECOSTRESSGriddedGranule:
        if overwrite is None:
            overwrite = cls._DEFAULT_OVERWRITE

        if output_filename is None and output_directory is None:
            output_directory = dirname(swath_granule.product_filename)
            output_filename = join(output_directory, cls.gridded_filename_base(swath_granule))

        if output_filename is None and output_directory is not None:
            output_filename = join(output_directory, cls.gridded_filename_base(swath_granule))

        if exists(output_filename):
            if overwrite:
                logger.info(f"existing file will be overwritten: {cl.file(output_filename)}")
            else:
                logger.info(f"file already exists: {cl.file(output_filename)}")
                return cls(output_filename)

        source_geometry = swath_granule.geometry

        if projection_system is None:
            projection_system = cls._DEFAULT_PROJECTION_SYSTEM

        if projection_system == "local_UTM":
            if cell_size is None:
                cell_size = cls._DEFAULT_UTM_CELL_SIZE

            target_geometry = source_geometry.UTM(cell_size)

        elif projection_system == "global_geographic":
            if cell_size is None:
                cell_size = cls._DEFAULT_GEOGRAPHIC_CELL_SIZE

            target_geometry = source_geometry.geographic(cell_size)

        else:
            raise ValueError(f"unrecognized projection system: {projection_system}")

        if grid_name is None:
            grid_name = cls._GRID_NAME

        if compression is None:
            compression = cls._DEFAULT_HDF5_COMPRESSION

        if apply_scale is None:
            apply_scale = cls._DEFAULT_GRIDDED_PRODUCT_APPLY_SCALE

        if apply_cloud is None:
            apply_cloud = cls._DEFAULT_GRIDDED_PRODUCT_APPLY_CLOUD

        logger.info("processing ECOSTRESS swath product to gridded product")
        logger.info(f"swath input: {cl.file(swath_granule.product_filename)}")
        logger.info(f"gridded output: {cl.file(output_filename)}")

        makedirs(output_directory, exist_ok=True)

        output_filename_partial = f"{output_filename}.partial"

        if exists(output_filename_partial):
            logger.warning(f"removing previous partial file: {output_filename_partial}")

        with HDFEOS5(output_filename_partial, "w") as output_file:
            output_file.write_geometry(target_geometry, grid_name)
            logger.info("copying ECOSTRESS metadata")
            swath_granule.copy_metadata(output_file, output_product_metadata_group_name=cls._PRODUCT_METADATA_GROUP)
            kdtree_timer = Timer()

            logger.info(
                f"building K-D tree for L1B geolocation: {cl.file(swath_granule.L1B_GEO_filename)}"
            )

            kd_tree = source_geometry.kd_tree(target_geometry)
            logger.info(f"finished building K-D tree ({cl.time(kdtree_timer)})")
            file_timer = Timer()
            logger.info(f"creating HDF5 file: {cl.file(output_filename)}")
            output_group = output_file.data_group(grid_name)

            for dataset_name in swath_granule.variables:
                logger.info(f"ingesting dataset: {dataset_name}")
                dataset_timer = Timer()

                source_data = swath_granule.variable(
                    dataset_name,
                    apply_scale=apply_scale,
                    apply_cloud=apply_cloud,
                    gridded=True,
                    kd_tree=kd_tree
                )

                source_metadata = swath_granule.metadata(dataset_name)
                _Fillvalue = source_metadata["nodata"]
                long_name = source_metadata["name"]
                units = source_metadata["units"]
                kwargs = {}

                if apply_scale:
                    add_offset = 0
                    scale_factor = 1

                    if scaleoffset is not None:
                        logger.info(f"HDF5 scaleoffset: {scaleoffset}")
                        kwargs["scaleoffset"] = scaleoffset
                else:
                    add_offset = source_metadata["offset"]
                    scale_factor = source_metadata["scale"]

                if str(source_data.dtype).startswith("float") and scaleoffset is not None:
                    logger.info(f"HDF5 scaleoffset: {scaleoffset}")
                    kwargs["scaleoffset"] = scaleoffset

                logger.info(f"finished ingesting dataset: {cl.name(dataset_name)} ({cl.time(dataset_timer)})")
                write_timer = Timer()
                output_dataset_name = cls.translate_dataset_name(dataset_name)
                logger.info(f"creating dataset: {output_dataset_name} ({source_data.dtype})")

                dataset = output_file.write(
                    name=output_dataset_name,
                    data=source_data,
                    compression=compression,
                    fillvalue=_Fillvalue,
                    group=output_group,
                    **kwargs
                )

                logger.info(f"{dataset_name}._Fillvalue: {_Fillvalue}")
                dataset.attrs.create("_Fillvalue", _Fillvalue, dtype="float64")
                logger.info(f"{dataset_name}.add_offset: {add_offset}")
                dataset.attrs.create("add_offset", add_offset, dtype="float64")
                logger.info(f"{dataset_name}.scale_factor: {scale_factor}")
                dataset.attrs.create("scale_factor", scale_factor, dtype="float64")
                logger.info(f"{dataset_name}.long_name: {long_name}")
                dataset.attrs.create("long_name", long_name, dtype=f"|S{len(long_name)}")
                logger.info(f"{dataset_name}.units: {units}")
                dataset.attrs.create("units", units, dtype=f"|S{len(units)}")
                logger.info(f"finished writing dataset: {cl.name(dataset_name)} ({cl.time(write_timer)})")

        shutil.move(output_filename_partial, output_filename)

        logger.info(f"completed HDF5 file: {cl.file(output_filename)} ({cl.time(file_timer)})")
        granule = cls(output_filename)

        return granule

    def tile_granule_name(self, tile: str) -> str:
        granule_name = self.granule_name
        granule_name = granule_name.replace(self.level, self.tiled_level)
        parts = granule_name.split("_")
        before_tile = "_".join(parts[:-3])
        after_tile = "_".join(parts[-3:])
        granule_name = f"{before_tile}_{tile}_{after_tile}"

        return granule_name

    @property
    def tile_grids(self) -> gpd.GeoDataFrame:
        return self.sentinel.tile_grids(target_geometry=self.boundary)

    @property
    def tile_footprints(self) -> gpd.GeoDataFrame:
        return self.sentinel.tile_footprints(target_geometry=self.boundary)

    @property
    def tiled_output_variables(self) -> List[str]:
        return self.variables

    @property
    def primary_variable(self) -> str:
        if self._PRIMARY_VARIABLE is None:
            raise ValueError(f"no primary variable for {self.__class__.__name__}")

        return self._PRIMARY_VARIABLE

    @property
    def TileClass(self):
        return self._TILE_CLASS

    def to_tile(
            self,
            tile: str,
            tile_granule_directory: str = None,
            granule_name: str = None,
            geometry: RasterGeometry = None,
            compression: str = None,
            variables: List[str] = None,
            overwrite: bool = False,
            skip_blank: bool = True) -> ECOSTRESSTiledGranule:
        if geometry is None:
            geometry = sentinel_tile_grid.grid(tile, cell_size=70)

        if not self.geometry.intersects(geometry):
            raise ValueError(
                f"tile geometry does not intersect scene geometry\ntile:\n{geometry}\nscene:\n{self.geometry}")

        TileClass = self.TileClass

        tile_granule = TileClass.from_scene(
            gridded_granule=self,
            tile=tile,
            tile_granule_directory=tile_granule_directory,
            tile_granule_name=granule_name,
            geometry=geometry,
            compression=compression,
            variables=variables,
            overwrite=overwrite,
            skip_blank=skip_blank
        )

        return tile_granule

    def to_tiles(
            self,
            output_directory: str = None,
            compression: str = None,
            level: int = None,
            overwrite: bool = False,
            tiles: List[str] = None,
            variables: List[str] = None,
            write_zip: bool = True,
            write_browse: bool = True,
            clear_directories: bool = True,
            skip_blank: bool = True) -> List[ECOSTRESSTiledGranule]:
        if output_directory is None:
            output_directory = os.getcwd()

        if compression is None:
            compression = self._DEFAULT_GEOTIFF_COMPRESSION

        if level is None:
            level = self._DEFAULT_COMPRESSION_LEVEL

        if isinstance(tiles, str):
            tiles = [tiles]

        tiled_scene_directory = output_directory
        # tiled_scene_directory = join(directory, self.tiled_scene_name)
        logger.info(f"tiled scene directory: {cl.dir(tiled_scene_directory)}")
        makedirs(tiled_scene_directory, exist_ok=True)
        tile_df = self.sentinel.tile_grids(target_geometry=self.boundary, eliminate_redundancy=True)
        tile_names = list(sorted(set(tile_df.tile)))

        logger.info(
            f"matching Sentinel tiles for orbit {self.orbit:05d} "
            f"scene {self.scene:03d}: {', '.join(tile_names)}"
        )

        tile_df["name"] = tile_df.tile.apply(lambda tile: self.tile_granule_name(tile))
        tile_df = tile_df[["tile", "name", "grid"]]
        tile_df.sort_values(by="name", inplace=True)
        tile_df["geometry"] = tile_df["grid"]
        tile_df.reset_index(drop=True, inplace=True)
        tile_df = tile_df[["tile", "name", "geometry"]]

        if tiles is not None:
            tile_df = tile_df[tile_df.tile.apply(lambda tile: tile in tiles)]

        output_tiles = []

        for i, (tile_name, granule_name, geometry) in tile_df.iterrows():
            logger.info(
                f"started processing tile {tile_name} ({i + 1} / {len(tile_df)}) "
                f"level: {self.tiled_level} "
                f"granule: {granule_name}"
            )

            tile_start_time = process_time()
            tile_granule_directory = join(tiled_scene_directory, granule_name)
            zip_filename = f"{tile_granule_directory}.zip"
            browse_filename = f"{tile_granule_directory}.png"

            output_filenames = []

            if write_zip:
                output_filenames.append(zip_filename)

            if write_browse:
                output_filenames.append(browse_filename)

            if not overwrite:
                if all([exists(filename) for filename in output_filenames]):
                    logger.info(f"product files already exist:")

                    for filename in sorted(output_filenames):
                        logger.info(f"* {filename}")

                    continue

            logger.info(f"creating tile granule directory: {cl.dir(tile_granule_directory)}")
            makedirs(tile_granule_directory, exist_ok=True)

            if not exists(tile_granule_directory):
                raise IOError(f"unable to create tile granule directory: {tile_granule_directory}")

            try:
                tile = self.to_tile(
                    tile=tile_name,
                    tile_granule_directory=tile_granule_directory,
                    granule_name=granule_name,
                    geometry=geometry,
                    compression=compression,
                    variables=variables,
                    overwrite=overwrite,
                    skip_blank=skip_blank
                )
            except BlankOutput as e:
                logger.warning(e)

                if clear_directories:
                    logger.info(f"removing tile granule directory: {tile_granule_directory}")
                    shutil.rmtree(tile_granule_directory, ignore_errors=True)

                continue
            except Exception as e:
                raise(e)

            if write_zip:
                tile.write_zip(zip_filename)

            if write_browse:
                tile.write_browse_image(browse_filename)

            if clear_directories:
                logger.info(f"removing tile granule directory: {tile_granule_directory}")
                shutil.rmtree(tile_granule_directory, ignore_errors=True)

            output_tiles.append(tile)
            tile_end_time = process_time()
            tile_duration = tile_end_time - tile_start_time

            logger.info(
                f"finished processing target {tile_name} ({i + 1} / {len(output_tiles)}) "
                f"level: {self.tiled_level} "
                f"granule: {granule_name} "
                f"({tile_duration:0.2f}s)"
            )

        return output_tiles

    def inspect(self):
        with HDFEOS5(self.product_filename, "r") as file:
            file.inspect()


