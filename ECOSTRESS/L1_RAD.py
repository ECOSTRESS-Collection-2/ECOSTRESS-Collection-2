from __future__ import annotations

import json
import logging
import os
import shutil
from datetime import datetime
from glob import glob
from os import makedirs
from os.path import exists, abspath, expanduser, dirname, join, basename, splitext
from typing import Union, List, Any

import numpy as np
import rasters as rt
import colored_logging as cl
import he5py
import rasters
from ECOSTRESS.L2_CLOUD import ECOv001L2CLOUD
from ECOSTRESS.L2_LSTE import L2LSTE, ECOv002L2LSTE
from ECOSTRESS.exit_codes import InputFilesInaccessible
from HDFEOS5.HDFEOS5 import HDFEOS5, h5py_copy
from rasters import Raster, KDTree, RasterGeometry, MultiRaster, RasterGrid
from timer import Timer
from .ECOSTRESS_granule import ECOSTRESSGranule
from .ECOSTRESS_gridded_tiled_granule import ECOSTRESSGriddedGranule, ECOSTRESSTiledGranule
from .ECOSTRESS_swath_granule import ECOSTRESSSwathGranule

__author__ = "Gregory Halverson"

from .XML_metadata import write_XML_metadata
from .scan_resampling import generate_scan_kd_trees, resample_scan_by_scan, clip_tails

logger = logging.getLogger(__name__)

PRIMARY_VARIABLE = "false_color"
PRODUCT_METADATA_GROUP = "L1B_RADMetadata"
GRID_NAME = "ECO_L1CG_RAD_70m"

L1CG_RAD_SHORT_NAME = "ECO_L1CG_RAD"
L1CG_RAD_LONG_NAME = "ECOSTRESS Gridded Top of Atmosphere Calibrated Radiance Instantaneous L1CG Global 70 m"

L1CT_RAD_SHORT_NAME = "ECO_L1CT_RAD"
L1CT_RAD_LONG_NAME = "ECOSTRESS Tiled Top of Atmosphere Calibrated Radiance Instantaneous L1CT Global 70 m"

TILED_OUTPUT_VARIABLES = [
    "cloud",
    "water",
    "radiance_1",
    "data_quality_1",
    "radiance_2",
    "data_quality_2",
    "radiance_3",
    "data_quality_3",
    "radiance_4",
    "data_quality_4",
    "radiance_5",
    "data_quality_5"
]


class L1RADGranule(ECOSTRESSGranule):
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE
    _PRODUCT_METADATA_GROUP = PRODUCT_METADATA_GROUP

    def __init__(self, product_filename: str):
        ECOSTRESSGranule.__init__(self, product_filename=product_filename)

    def radiance(self, band: int, apply_cloud: bool = False) -> Raster:
        if band not in range(1, 6):
            raise ValueError(f"invalid radiance band: {band}")

        radiance = self.variable(f"radiance_{band}", apply_cloud=apply_cloud)

        return radiance

    @property
    def false_color(self) -> MultiRaster:
        band_5 = self.radiance(5, apply_cloud=False)
        band_4 = self.radiance(4, apply_cloud=False)
        band_2 = self.radiance(2, apply_cloud=False)
        false_color = MultiRaster.stack([band_5, band_4, band_2])

        return false_color

    def variable(
            self,
            variable_name: str,
            apply_scale: bool = True,
            apply_cloud: bool = True,
            geometry: RasterGeometry = None,
            fill_value: Any = None,
            kd_tree: KDTree = None,
            **kwargs):
        if self.data_group_name is None:
            raise ValueError("no data group name")

        if variable_name == "water":
            data = self.water.astype(np.uint8)
            fill_value = 0
        elif variable_name == "cloud":
            data = self.cloud.astype(np.uint8)
            fill_value = 0
        elif variable_name == "false_color":
            data = self.false_color.astype(np.float32)
        else:
            if "quality" in variable_name:
                apply_scale = False
                apply_cloud = False

            data = self.data(
                dataset_name=f"{self.data_group_name}/{variable_name}",
                apply_scale=apply_scale,
                apply_cloud=apply_cloud
            )

        if geometry is not None:
            gridded_data = data.resample(target_geometry=geometry, nodata=fill_value, kd_tree=kd_tree)
            data = gridded_data.astype(data.dtype)

        return data


class L1BRAD(ECOSTRESSSwathGranule, L1RADGranule):
    logger = logging.getLogger(__name__)

    _DATA_GROUP = "Radiance"
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE
    _PRODUCT_METADATA_GROUP = PRODUCT_METADATA_GROUP

    def __init__(
            self,
            L1B_RAD_filename: str,
            L2_CLOUD_filename: str = None,
            L2_LSTE_filename: str = None,
            L1B_GEO_filename: str = None,
            **kwargs):
        L1B_RAD_filename = abspath(expanduser(L1B_RAD_filename))

        orbit = int(splitext(basename(L1B_RAD_filename))[0].split("_")[-5])
        scene = int(splitext(basename(L1B_RAD_filename))[0].split("_")[-4])

        if L2_CLOUD_filename is None:
            directory = abspath(expanduser(dirname(L1B_RAD_filename)))
            pattern = join(directory, f"*_L2_CLOUD_{orbit:05d}_{scene:03d}_*.h5")
            logger.info(f"searching for L2 CLOUD: {cl.val(pattern)}")
            candidates = sorted(glob(pattern))

            if len(candidates) == 0:
                raise ValueError("no L2 CLOUD filename given or found")

            L2_CLOUD_filename = candidates[-1]

            logger.info(f"found L2 CLOUD file: {cl.file(L2_CLOUD_filename)}")

        if L1B_GEO_filename is None:
            directory = abspath(expanduser(dirname(L1B_RAD_filename)))
            pattern = join(directory, f"*_L1B_GEO_{orbit:05d}_{scene:03d}_*.h5")
            logger.info(f"searching for L1B GEO: {cl.val(pattern)}")
            candidates = sorted(glob(pattern))

            if len(candidates) == 0:
                raise ValueError("no L1B GEO filename given or found")

            L1B_GEO_filename = candidates[-1]

            logger.info(f"found L1B GEO file: {cl.file(L1B_GEO_filename)}")

        L2_CLOUD_filename = abspath(expanduser(L2_CLOUD_filename))
        L1B_GEO_filename = abspath(expanduser(L1B_GEO_filename))

        ECOSTRESSSwathGranule.__init__(
            self,
            product_filename=L1B_RAD_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            data_group_name=self._DATA_GROUP,
            product_metadata_group_name=self._PRODUCT_METADATA_GROUP,
            **kwargs
        )

        L1RADGranule.__init__(self, product_filename=L1B_RAD_filename)

        if not exists(L1B_RAD_filename):
            raise InputFilesInaccessible(f"L1 RAD file does not exist: {L1B_RAD_filename}")

        self.L2_LSTE_filename = None
        self.L2_LSTE = None

        if L2_LSTE_filename is not None:
            self.L2_LSTE_filename = L2_LSTE_filename
            self.L2_LSTE = L2LSTE.open(
                L2_LSTE_filename=L2_LSTE_filename,
                L2_CLOUD_filename=L2_CLOUD_filename,
                L1B_GEO_filename=L1B_GEO_filename
            )

    def __repr__(self):
        display_dict = {
            "L1B RAD": self.product_filename,
            "L2 CLOUD:": self.L2_CLOUD_filename,
            "L1B GEO:": self.L1B_GEO_filename
        }

        display_string = f"L1B RAD Granule\n{json.dumps(display_dict, indent=2)}"

        return display_string

    @property
    def water(self) -> Raster:
        if isinstance(self.L2_LSTE, ECOv002L2LSTE):
            return self.L2_LSTE.water
        elif isinstance(self.L2_CLOUD, ECOv001L2CLOUD):
            return self.L2_CLOUD.water
        else:
            raise ValueError(
                f"unable to retrieve water mask with L2 LSTE type {type(self.L2_LSTE)} and L2 CLOUD type {type(self.L2_CLOUD)}")

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

        nodata = np.nan

        if variable_name == "view_zenith":
            data = self.L1B_GEO.view_zenith.astype(np.float32)
        elif variable_name == "height":
            data = self.L1B_GEO.height.astype(np.float32)
        elif variable_name == "water":
            data = self.water.astype(np.uint8)
        elif variable_name == "cloud":
            data = self.cloud.astype(np.uint8)
        elif variable_name == "false_color":
            data = self.false_color.astype(np.float32)
        else:
            data = self.data(dataset_name=f"{self.data_group_name}/{variable_name}", apply_scale=apply_scale, apply_cloud=apply_cloud)

            if apply_scale:
                logger.info("filtering values less than -9990")
                data = rt.where(data <= -9990, np.nan, data)

            logger.info(f"{variable_name} min: {np.nanmin(data)} mean: {np.nanmean(data)} max: {np.nanmax(data)}")

        dtype = data.dtype

        if "uint8" in str(dtype):
            # nodata = 255
            nodata = 0
        elif "uint16" in str(dtype):
            nodata = 65535

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
                gridded_data = clip_tails(data).to_geometry(geometry, nodata=nodata, kd_tree=kd_tree)
                logger.info(f"finished resampling {variable_name} with clipped-tails strategy ({timer})")
            else:
                raise ValueError(f"unrecognized overlap strategy: {overlap_strategy}")

            data = gridded_data.astype(dtype)
            data.nodata = nodata

        if "cloud" in variable_name or "water" in variable_name:
            data = rt.where(data == 1, 1, 0)
            data = data.astype(bool)

        return data

    def metadata(self, variable_name: str) -> dict:
        if variable_name in ["water", "cloud", "view_zenith"]:
            return {}
        else:
            return super(L1BRAD, self).metadata(variable_name=variable_name)

    @property
    def variables(self) -> List[str]:
        return ["water", "cloud"] + super().variables


class L1CTRAD(ECOSTRESSTiledGranule, L1RADGranule):
    _PRODUCT_NAME = "L1CT_RAD"
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
            containing_directory: str = None,
            **kwargs):
        L1RADGranule.__init__(self, product_filename=product_location)
        ECOSTRESSTiledGranule.__init__(
            self,
            product_location=product_location,
            orbit=orbit,
            scene=scene,
            tile=tile,
            time_UTC=time_UTC,
            build=build,
            process_count=process_count,
            containing_directory=containing_directory,
            **kwargs
        )

        self._water = None
        self._cloud = None

    @property
    def filename(self) -> str:
        return self.product_location

    @property
    def geometry(self) -> RasterGrid:
        if self._geometry is not None:
            return self._geometry

        URI = self.layer_URI("radiance_5")
        grid = RasterGrid.open(URI)

        return grid

    @property
    def variables(self) -> List[str]:
        return [
            "radiance_1",
            "radiance_2",
            "radiance_3",
            "radiance_4",
            "radiance_5",
            "data_quality_1",
            "data_quality_2",
            "data_quality_3",
            "data_quality_4",
            "data_quality_5"
            "water",
            "cloud"
        ]

    @property
    def water(self) -> Raster:
        if self._water is None:
            self._water = self.variable(variable="water")

        return self._water

    @property
    def cloud(self) -> Raster:
        if self._cloud is None:
            self._cloud = self.variable(variable="cloud")

        return self._cloud

    def variable(self, variable: str, geometry: RasterGeometry = None, **kwargs) -> Raster:
        if variable == "false_color":
            return self.false_color

        return ECOSTRESSTiledGranule.variable(self, variable=variable, geometry=geometry)

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
            skip_blank: bool = False) -> L1CTRAD:
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
        standard_metadata = metadata[cls._STANDARD_METADATA_GROUP]

        DATA_FORMAT_TYPE = "COG"
        logger.info(f"DataFormatType: {cl.val(DATA_FORMAT_TYPE)}")
        standard_metadata["DataFormatType"] = DATA_FORMAT_TYPE

        standard_metadata.pop("HDFVersionID")

        PROCESSING_LEVEL = "L1CT"
        logger.info(f"ProcessingLevelID: {cl.val(PROCESSING_LEVEL)}")
        standard_metadata["ProcessingLevelID"] = PROCESSING_LEVEL

        PROCESSING_LEVEL_DESCRIPTION = "Level 1 Tiled Top of Atmosphere Calibrated Radiance"
        logger.info(f"ProcessingLevelDescription: {cl.val(PROCESSING_LEVEL_DESCRIPTION)}")
        standard_metadata["ProcessingLevelDescription"] = PROCESSING_LEVEL_DESCRIPTION

        SIS_VERSION = "Preliminary"
        logger.info(f"SISVersion: {cl.val(SIS_VERSION)}")
        standard_metadata["SISVersion"] = SIS_VERSION

        short_name = L1CT_RAD_SHORT_NAME
        logger.info(f"L1CT RAD short name: {cl.val(short_name)}")
        standard_metadata["ShortName"] = short_name

        long_name = L1CT_RAD_LONG_NAME
        logger.info(f"L1CT RAD long name: {cl.val(long_name)}")
        standard_metadata["LongName"] = long_name

        logger.info(f"RegionID: {cl.place(tile)}")
        standard_metadata["RegionID"] = tile

        cell_width = geometry.cell_width
        rows, cols = geometry.shape
        bbox = geometry.bbox.latlon
        x_min, y_min, x_max, y_max = bbox
        CRS = geometry.proj4

        logger.info(f"ImageLineSpacing: {cl.val(cell_width)}")
        standard_metadata["ImageLineSpacing"] = cell_width
        logger.info(f"ImagePixelSpacing: {cl.val(cell_width)}")
        standard_metadata["ImagePixelSpacing"] = cell_width
        logger.info(f"ImageLines: {cl.val(rows)}")
        standard_metadata["ImageLines"] = rows
        logger.info(f"ImagePixels: {cl.val(cols)}")
        standard_metadata["ImagePixels"] = cols
        logger.info(f"CRS: {cl.val(CRS)}")
        standard_metadata["CRS"] = CRS
        logger.info(f"EastBoundingCoordinate: {cl.place(x_max)}")
        standard_metadata["EastBoundingCoordinate"] = x_max
        logger.info(f"SouthBoundingCoordinate: {cl.place(y_min)}")
        standard_metadata["SouthBoundingCoordinate"] = y_min
        logger.info(f"WestBoundingCoordinate: {cl.place(x_min)}")
        standard_metadata["WestBoundingCoordinate"] = x_min
        logger.info(f"NorthBoundingCoordinate: {cl.place(y_max)}")
        standard_metadata["NorthBoundingCoordinate"] = y_max

        product_metadata = metadata[cls._PRODUCT_METADATA_GROUP]

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
                apply_cloud=False,
                geometry=geometry
            )

            granule.add_layer(variable, image)

            logger.info(
                f"finished processing variable {variable} ({j + 1} / {len(output_variables)}) "
                f"for granule: {tile_granule_name} "
                f"({cl.time(timer)})"
            )

        cloud = granule.cloud
        cloud_percent = np.count_nonzero(cloud) / cloud.size * 100

        logger.info(f"QAPercentCloudCover: {cl.val(cloud_percent)}")
        product_metadata["QAPercentCloudCover"] = cloud_percent

        metadata[cls._STANDARD_METADATA_GROUP] = standard_metadata
        metadata[cls._PRODUCT_METADATA_GROUP] = product_metadata
        granule.write_metadata(metadata)

        false_color = granule.false_color
        granule.add_layer("false_color", false_color, include_COG=False, include_geojpeg=True)

        return granule


class L1CGRAD(ECOSTRESSGriddedGranule, L1RADGranule):
    _TILE_CLASS = L1CTRAD
    _GRID_NAME = GRID_NAME
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE
    _PRODUCT_METADATA_GROUP = "ProductMetadata"

    TILED_OUTPUT_VARIABLES = TILED_OUTPUT_VARIABLES

    def __init__(self, L1CG_RAD_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L1CG_RAD_filename)
        L1RADGranule.__init__(self, product_filename=L1CG_RAD_filename)

    # @property
    # def tiled_output_variables(self):
    #     tiled_output_variables = ["cloud", "water"]
    #
    #     for band in self.available_bands:
    #         tiled_output_variables.append(f"radiance_{band}")
    #         tiled_output_variables.append(f"data_quality_{band}")
    #
    #     return tiled_output_variables

    @property
    def cloud(self) -> Raster:
        return self.variable(
            variable_name="cloud",
            apply_scale=False,
            apply_cloud=False
        )

    @property
    def water(self) -> Raster:
        return self.variable(
            variable_name="water",
            apply_scale=False,
            apply_cloud=False
        )

    def variable(
            self,
            variable_name: str,
            apply_scale: bool = True,
            apply_cloud: bool = True,
            geometry: RasterGeometry = None,
            **kwargs) -> Raster:

        if variable_name == "false_color":
            return self.false_color

        if variable_name in ["cloud", "water"] or variable_name.startswith("data_quality"):
            apply_scale = False
            apply_cloud = False

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

        if apply_cloud:
            cloud = self.variable(
                variable_name="cloud",
                apply_scale=False,
                apply_cloud=False,
                geometry=geometry
            )

            image = rasters.where(cloud, np.nan, image)

        return image

    @classmethod
    def from_swath(
            cls,
            swath_granule: L1BRAD,
            output_filename: str = None,
            output_directory: str = None,
            projection_system: str = None,
            cell_size_degrees: float = 0.0007,
            cell_size_meters: float = 70,
            search_radius_meters: float = 100,
            grid_name: str = None,
            compression: str = None,
            scaleoffset: int = None,
            gridded_geometry: RasterGrid = None,
            overlap_strategy: str = "checkerboard",
            kd_tree: KDTree = None,
            scan_kd_trees: List[KDTree] = None,
            kd_tree_path: str = None,
            variables: List[str] = None,
            PGE_name: str = None,
            PGE_version: str = None,
            build: str = None,
            input_filenames: List[str] = None,
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

        swath_geometry = swath_granule.geometry

        if projection_system is None:
            projection_system = cls._DEFAULT_PROJECTION_SYSTEM

        if projection_system == "local_UTM":
            logger.info("using local UTM projection system for gridded products")

            if gridded_geometry is None:
                gridded_geometry = swath_geometry.UTM(cell_size_meters)

            cell_size = cell_size_meters
        elif projection_system == "global_geographic":
            logger.info("using global geographic projection system for gridded products")

            if gridded_geometry is None:
                gridded_geometry = swath_geometry.geographic(cell_size_degrees)

            cell_size = cell_size_degrees
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

        if kd_tree_path is not None and not exists(kd_tree_path):
            logger.warning(f"K-D tree file not found: {kd_tree_path}")

        # TODO consolidate redundancy in preparing K-D trees between L1_RAD.py, L2_LSTE.py, L2G_CLOUD.py, L1_L1_RAD_LSTE.py

        if overlap_strategy == "checkerboard":
            logger.info("using checkerboard overlap strategy")

            # treating K-D tree path as filename for whole-scene processing
            kd_tree_filename = kd_tree_path

            if kd_tree_filename is not None and exists(kd_tree_filename):
                logger.info(f"started loading checkerboard K-D tree: {kd_tree_filename}")
                timer = Timer()
                kd_tree = KDTree.load(kd_tree_filename)
                logger.info(f"finished loading checkerboard K-D tree: {kd_tree_filename} ({timer})")
            else:
                logger.info("started building checkerboard K-D tree")
                timer = Timer()

                kd_tree = KDTree(
                    source_geometry=swath_geometry,
                    target_geometry=gridded_geometry,
                    radius_of_influence=search_radius_meters
                )

                logger.info(f"finished building checkerboard K-D tree ({timer})")

                if kd_tree_filename is not None:
                    logger.info(f"started saving checkerboard K-D tree: {kd_tree_filename}")
                    timer = Timer()
                    kd_tree.save(kd_tree_filename)
                    logger.info(f"finished saving checkerboard K-D tree ({timer}): {kd_tree_filename}")

        elif overlap_strategy == "scan_by_scan":
            logger.info("using scan-by-scan overlap strategy")

            # treating the K-D tree path as a directory for scan-by-scan approach
            kd_tree_directory = kd_tree_path

            if scan_kd_trees is None:
                if kd_tree_directory is not None and exists(kd_tree_directory):
                    logger.info("started loading scan-by-scan K-D trees")
                    timer = Timer()

                    scan_kd_trees = [KDTree.load(filename) for filename in sorted(glob(join(kd_tree_directory, "*.kdtree")))]

                    logger.info(f"finished loading scan-by-scan K-D trees ({cl.time(timer)})")
                else:
                    logger.info("started building scan-by-scan K-D trees")
                    timer = Timer()

                    scan_kd_trees = generate_scan_kd_trees(
                        swath_geometry=swath_geometry,
                        cell_size_degrees=cell_size_degrees
                    )

                    logger.info(f"finished building scan-by-scan K-D trees ({cl.time(timer)})")

                    if kd_tree_directory is not None:
                        makedirs(kd_tree_directory, exist_ok=True)

                        for i, kd_tree in enumerate(scan_kd_trees):
                            # formatting K-D tree filenames with two-digit leading-zero enumeration for sorting in scan-by-scan approach
                            kd_tree_filename = join(kd_tree_directory, f"{i:02d}.kdtree")
                            kd_tree.save(kd_tree_filename)

        elif overlap_strategy == "remove_105_128":
            logger.info("using clipped-tails overlap strategy")

            kd_tree_filename = kd_tree_path

            logger.info(f"original swath shape: {swath_geometry.shape}")
            swath_geometry = clip_tails(swath_geometry)
            logger.info(f"clipped swath shape: {swath_geometry.shape}")

            if kd_tree_filename is not None and exists(kd_tree_filename):
                logger.info(f"started loading checkerboard K-D tree: {cl.file(kd_tree_filename)}")
                timer = Timer()
                kd_tree = KDTree.load(kd_tree_filename)
                logger.info(f"finished loading checkerboard K-D tree: {cl.file(kd_tree_filename)} ({cl.time(timer)})")
            else:
                logger.info("started building clipped-tails K-D tree")
                timer = Timer()

                kd_tree = KDTree(
                    source_geometry=swath_geometry,
                    target_geometry=gridded_geometry,
                    radius_of_influence=search_radius_meters
                )

                logger.info(f"finished building clipped-tails K-D tree ({cl.time(timer)})")

                if kd_tree_filename is not None:
                    logger.info(f"started saving clipped-tails K-D tree: {cl.file(kd_tree_filename)}")
                    timer = Timer()
                    kd_tree.save(kd_tree_filename)
                    logger.info(
                        f"finished saving clipped-tails K-D tree ({cl.time(timer)}): {cl.file(kd_tree_filename)}")
        else:
            raise ValueError(f"unrecognized overlap strategy: {overlap_strategy}")

        makedirs(dirname(output_filename), exist_ok=True)
        dataset_metadata = {}
        output_filename_partial = f"{output_filename}.partial"
        output_filename_rewritten = f"{output_filename}.rewritten"

        if exists(output_filename_partial):
            os.remove(output_filename_partial)

        if variables is None:
            variables = swath_granule.variables

        with he5py.File(output_filename_partial, "w") as output_file:
            file_timer = Timer()
            logger.info(f"creating HDF-EOS5 file: {cl.file(output_filename_partial)}")

            output_grid = output_file.create_grid(grid_name=grid_name, geometry=gridded_geometry)

            for dataset_name in variables:
                logger.info(f"ingesting dataset: {dataset_name}")
                dataset_timer = Timer()

                source_data = swath_granule.variable(
                    dataset_name,
                    apply_scale=True,
                    apply_cloud=False,
                    geometry=gridded_geometry,
                    kd_tree=kd_tree,
                    scan_kd_trees=scan_kd_trees,
                    overlap_strategy=overlap_strategy
                )

                source_metadata = swath_granule.metadata(dataset_name)

                if "nodata" in source_metadata:
                    _Fillvalue = source_metadata["nodata"]
                else:
                    _Fillvalue = None

                if "name" in source_metadata:
                    long_name = source_metadata["name"]
                else:
                    long_name = dataset_name

                if "units" in source_metadata:
                    units = source_metadata["units"]
                else:
                    units = "unitless"

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
                # logger.info(f"creating dataset: {output_dataset_name} ({source_data.dtype})")

                if source_data.dtype in (np.float32, np.float64):
                    source_data = source_data.astype(np.float32)
                    _Fillvalue = np.nan

                    logger.info(f"creating dataset: {cl.name(output_dataset_name)} ({cl.val(source_data.dtype)})")

                    output_grid.write_float(
                        field_name=output_dataset_name,
                        image=source_data
                    )
                elif source_data.dtype in (np.uint8, np.int8, bool):
                    source_data = source_data.astype(np.uint8)
                    logger.info(f"creating dataset: {cl.name(output_dataset_name)} ({cl.val(source_data.dtype)})")

                    if _Fillvalue is None:
                        _Fillvalue = cls.UINT8_FILL

                    output_grid.write_uint8(
                        field_name=output_dataset_name,
                        image=source_data
                    )
                elif source_data.dtype in (np.uint16, np.int16):
                    source_data = source_data.astype(np.uint16)
                    logger.info(f"creating dataset: {cl.name(output_dataset_name)} ({cl.val(source_data.dtype)})")

                    if _Fillvalue is None:
                        _Fillvalue = cls.UINT16_FILL

                    output_grid.write_float(
                        field_name=output_dataset_name,
                        image=source_data.astype(np.float32)
                    )
                else:
                    raise ValueError(f"unsupported source dtype: {source_data.dtype}")

                logger.info(f"finished writing dataset: {cl.name(dataset_name)} ({cl.time(write_timer)})")

                dataset_metadata[dataset_name] = {
                    "_Fillvalue": _Fillvalue,
                    "add_offset": add_offset,
                    "scale_factor": scale_factor,
                    "long_name": long_name,
                    "units": units
                }

        with HDFEOS5(output_filename_partial, "r+") as output_file:
            file_timer = Timer()
            logger.info(f"started copying ECOSTRESS metadata to HDF-EOS5 file: {cl.file(output_filename_partial)}")
            swath_granule.copy_metadata(output_file, output_product_metadata_group_name=cls._PRODUCT_METADATA_GROUP)
            standard_metadata = output_file[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{cls._STANDARD_METADATA_GROUP}"]

            HDF_version = output_file.libver[1]
            logger.info(f"HDFVersionID: {cl.val(HDF_version)}")
            standard_metadata["HDFVersionID"][...] = str(HDF_version)

            standard_metadata["AncillaryInputPointer"][...] = ""
            orbit = swath_granule.orbit
            scene = swath_granule.scene
            standard_metadata["RegionID"][...] = f"{orbit:05d}_{scene:03d}"

            short_name = L1CG_RAD_SHORT_NAME
            logger.info(f"L1CG RAD short name: {cl.val(short_name)}")
            standard_metadata["ShortName"][...] = short_name

            long_name = L1CG_RAD_LONG_NAME
            logger.info(f"L1CG RAD long name: {cl.val(long_name)}")
            standard_metadata["LongName"][...] = long_name

            if PGE_name is not None:
                logger.info(f"PGEName: {cl.name(PGE_name)}")
                standard_metadata["PGEName"][...] = str(PGE_name)

            if PGE_version is not None:
                logger.info(f"PGEVersion: {cl.name(PGE_version)}")
                standard_metadata["PGEVersion"][...] = str(PGE_version)

            boundary_WKT = swath_granule.boundary_WKT
            boundary_WKT_NAME = cls._BOUNDARY_WKT_NAME
            logger.info(f"{boundary_WKT_NAME}: {cl.place(boundary_WKT)}")
            standard_metadata.create_dataset(boundary_WKT_NAME, data=boundary_WKT)

            CRS = gridded_geometry.proj4
            logger.info(f"CRS: {cl.val(CRS)}")
            standard_metadata.create_dataset("CRS", data=CRS)

            bbox = gridded_geometry.bbox.latlon
            lon_min, lat_min, lon_max, lat_max = bbox

            logger.info(f"EastBoundingCoordinate: {cl.place(lon_max)}")
            standard_metadata["EastBoundingCoordinate"][...] = lon_max
            logger.info(f"SouthBoundingCoordinate: {cl.place(lat_min)}")
            standard_metadata["SouthBoundingCoordinate"][...] = lat_min
            logger.info(f"WestBoundingCoordinate: {cl.place(lon_min)}")
            standard_metadata["WestBoundingCoordinate"][...] = lon_min
            logger.info(f"NorthBoundingCoordinate: {cl.place(lat_max)}")
            standard_metadata["NorthBoundingCoordinate"][...] = lat_max

            DATA_FORMAT_TYPE = "HDF-EOS5"
            logger.info(f"DataFormatType: {cl.val(DATA_FORMAT_TYPE)}")
            standard_metadata["DataFormatType"][...] = str(DATA_FORMAT_TYPE)

            COLLECTION_LABEL = "ECOv002"
            logger.info(f"CollectionLabel: {cl.val(COLLECTION_LABEL)}")
            standard_metadata["CollectionLabel"][...] = str(COLLECTION_LABEL)

            if build is not None:
                logger.info(f"BuildID: {cl.val(build)}")
                standard_metadata["BuildID"][...] = str(build)

            rows, cols = gridded_geometry.shape

            logger.info(f"ImageLineSpacing: {cl.val(cell_size_degrees)}")
            standard_metadata["ImageLineSpacing"][...] = cell_size_degrees
            logger.info(f"ImagePixelSpacing: {cl.val(cell_size_degrees)}")
            standard_metadata["ImagePixelSpacing"][...] = cell_size_degrees
            logger.info(f"ImageLines: {cl.val(rows)}")
            standard_metadata["ImageLines"][...] = rows
            logger.info(f"ImagePixels: {cl.val(cols)}")
            standard_metadata["ImagePixels"][...] = cols

            standard_metadata["AncillaryInputPointer"][...] = ""

            if input_filenames is not None:
                input_pointer = ",".join([basename(filename) for filename in input_filenames])
                logger.info(f"InputPointer: {cl.place(input_pointer)}")
                standard_metadata["InputPointer"][...] = input_pointer

            local_granule_ID = basename(output_filename)
            logger.info(f"LocalGranuleID: {cl.file(local_granule_ID)}")
            standard_metadata["LocalGranuleID"][...] = str(local_granule_ID)

            PROCESSING_LEVEL = "L1CG"
            LONG_NAME = "Top of Atmosphere Calibrated Radiance"
            PROCESSING_LEVEL_DESCRIPTION = "Level 1 Gridded Top of Atmosphere Calibrated Radiance"
            logger.info(f"ProcessingLevelID: {cl.file(PROCESSING_LEVEL)}")
            standard_metadata["ProcessingLevelID"][...] = str(PROCESSING_LEVEL)
            standard_metadata["ProcessingLevelDescription"][...] = str(PROCESSING_LEVEL_DESCRIPTION)

            # product_metadata = output_file[f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{cls._PRODUCT_METADATA_GROUP}"]

            output_group = output_file.data_group(grid_name)

            for dataset_name in variables:
                logger.info(f"ingesting dataset: {dataset_name}")
                dataset_timer = Timer()
                dataset = output_group[dataset_name]

                _Fillvalue = dataset_metadata[dataset_name]["_Fillvalue"]
                add_offset = dataset_metadata[dataset_name]["add_offset"]
                scale_factor = dataset_metadata[dataset_name]["scale_factor"]
                long_name = dataset_metadata[dataset_name]["long_name"]
                units = dataset_metadata[dataset_name]["units"]

                if "QC" in dataset_name or "quality" in dataset_name:
                    dtype = "uint16"
                else:
                    dtype = str(dataset.dtype)

                logger.info(f"dataset dtype: {dtype}")

                if dtype == "uint8":
                    _Fillvalue = cls.UINT8_FILL
                elif dtype == "uint16":
                    _Fillvalue = cls.UINT16_FILL

                logger.info(f"{dataset_name}._Fillvalue ({dtype}): {_Fillvalue}")
                dataset.attrs.create("_Fillvalue", _Fillvalue, dtype=dtype)

                logger.info(f"{dataset_name}.add_offset: {add_offset}")
                dataset.attrs.create("add_offset", add_offset, dtype="float64")
                logger.info(f"{dataset_name}.scale_factor: {scale_factor}")
                dataset.attrs.create("scale_factor", scale_factor, dtype="float64")
                logger.info(f"{dataset_name}.long_name: {long_name}")
                dataset.attrs.create("long_name", long_name, dtype=f"|S{len(long_name)}")
                logger.info(f"{dataset_name}.units: {units}")
                dataset.attrs.create("units", units, dtype=f"|S{len(units)}")

                logger.info(f"finished ingesting dataset: {cl.name(dataset_name)} ({cl.time(dataset_timer)})")

            logger.info(
                f"finished copying ECOSTRESS metadata to HDF-EOS5 file: {cl.file(output_filename_partial)} ({cl.time(file_timer)})")

        logger.info(f"re-writing HDF-EOS5 file: {output_filename_partial} -> {output_filename_rewritten}")
        h5py_copy(output_filename_partial, output_filename_rewritten)
        logger.info(f"removing original HDF-EOS5 output: {output_filename_partial}")
        os.remove(output_filename_partial)
        logger.info(f"renaming completed output file: {output_filename_rewritten} -> {output_filename}")
        shutil.move(output_filename_rewritten, output_filename)
        logger.info(f"completed HDF-EOS5 generation: {cl.file(output_filename)}")
        granule = cls(output_filename)
        XML_metadata_filename = output_filename.replace(".h5", ".h5.met")
        logger.info(f"writing XML metadata file: {cl.file(XML_metadata_filename)}")
        write_XML_metadata(granule.standard_metadata, XML_metadata_filename)

        return granule
