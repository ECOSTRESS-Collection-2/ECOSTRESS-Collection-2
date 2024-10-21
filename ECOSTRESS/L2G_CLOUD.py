from __future__ import annotations

import logging
import os
import shutil
from datetime import datetime
from os.path import exists, basename, dirname, join
from typing import Union, List
from glob import glob
import PIL
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import Colormap
from os import makedirs
import colored_logging as cl
import he5py
import rasters as rt
from ECOSTRESS.ECOSTRESS_gridded_tiled_granule import ECOSTRESSGriddedGranule, ECOSTRESSTiledGranule
from ECOSTRESS.L2_CLOUD import L2CLOUD, L2CLOUDGranule
from ECOSTRESS.XML_metadata import write_XML_metadata
from ECOSTRESS.exit_codes import BlankOutput
from ECOSTRESS.scan_resampling import generate_scan_kd_trees, clip_tails
from ECOSTRESS_colors import CLOUD_CMAP
from HDFEOS5.HDFEOS5 import HDFEOS5, h5py_copy
from rasters import Raster, RasterGeometry, KDTree, RasterGrid
from timer import Timer

PRIMARY_VARIABLE = "Cloud_test_final"
PRODUCT_METADATA_GROUP = "L2 CLOUD Metadata"
GRID_NAME = "ECO_L2G_CLOUD_70m"

L2G_CLOUD_SHORT_NAME = "ECO_L2G_CLOUD"
L2G_CLOUD_LONG_NAME = "ECOSTRESS Gridded Cloud Mask Instantaneous L2 Global 70 m"

L2T_CLOUD_SHORT_NAME = "ECO_L2T_CLOUD"
L2T_CLOUD_LONG_NAME = "ECOSTRESS Tiled Cloud Mask Instantaneous L2 Global 70 m"

VARIABLE_NAMES = [
    "Cloud_test_1",
    "Cloud_test_2",
    "Cloud_test_final"
]

logger = logging.getLogger(__name__)


class L2TCLOUD(ECOSTRESSTiledGranule, L2CLOUD):
    _PRODUCT_NAME = "L2T_CLOUD"
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
        L2CLOUD.__init__(self, product_filename=product_location)
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
            skip_blank: bool = True) -> L2TCLOUD:
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

        PROCESSING_LEVEL = "L2T"
        logger.info(f"ProcessingLevelID: {cl.val(PROCESSING_LEVEL)}")
        standard_metadata["ProcessingLevelID"] = PROCESSING_LEVEL
        PROCESSING_LEVEL_DESCRIPTION = "Level 2 Tiled Cloud Mask"
        logger.info(f"ProcessingLevelDescription: {cl.val(PROCESSING_LEVEL_DESCRIPTION)}")
        standard_metadata["ProcessingLevelDescription"] = PROCESSING_LEVEL_DESCRIPTION
        SIS_VERSION = "Preliminary"
        logger.info(f"SISVersion: {cl.val(SIS_VERSION)}")
        standard_metadata["SISVersion"] = SIS_VERSION

        short_name = L2T_CLOUD_SHORT_NAME
        logger.info(f"L2T CLOUD short name: {cl.val(short_name)}")
        standard_metadata["ShortName"] = short_name

        long_name = L2T_CLOUD_LONG_NAME
        logger.info(f"L2T CLOUD long name: {cl.val(long_name)}")
        standard_metadata["LongName"] = long_name

        logger.info(f"RegionID: {cl.place(tile)}")
        standard_metadata["RegionID"] = tile

        cell_width = geometry.cell_width
        rows, cols = geometry.shape
        bbox = geometry.bbox
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

        product_metadata.pop("CloudMaxTemperature")
        product_metadata.pop("CloudMeanTemperature")
        product_metadata.pop("CloudMinTemperature")
        product_metadata.pop("CloudSDevTemperature")
        product_metadata.pop("Emis1GoodAvg")
        product_metadata.pop("Emis2GoodAvg")
        product_metadata.pop("Emis3GoodAvg")
        product_metadata.pop("Emis4GoodAvg")
        product_metadata.pop("Emis5GoodAvg")
        product_metadata.pop("LSTGoodAvg")

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

        cloud = granule.cloud
        cloud_percent = np.count_nonzero(cloud) / cloud.size * 100

        logger.info(f"QAPercentCloudCover: {cl.val(cloud_percent)}")
        product_metadata["QAPercentCloudCover"] = cloud_percent

        ST_K = granule.ST_K
        good_percent = np.count_nonzero(~np.isnan(ST_K)) / np.array(ST_K).size * 100
        logger.info(f"QAPercentGoodQuality: {cl.val(good_percent)}")
        product_metadata["QAPercentGoodQuality"] = good_percent

        metadata[cls._STANDARD_METADATA_GROUP] = standard_metadata
        metadata[cls._PRODUCT_METADATA_GROUP] = product_metadata
        granule.write_metadata(metadata)

        return granule


class L2GCLOUD(ECOSTRESSGriddedGranule, L2CLOUDGranule):
    _TILE_CLASS = L2TCLOUD
    _GRID_NAME = GRID_NAME
    _PRIMARY_VARIABLE = PRIMARY_VARIABLE
    _PRODUCT_METADATA_GROUP = "ProductMetadata"
    _GRANULE_PREVIEW_CMAP = CLOUD_CMAP

    VARIABLE_NAMES = VARIABLE_NAMES

    def __init__(self, L2G_CLOUD_filename: str):
        ECOSTRESSGriddedGranule.__init__(self, product_filename=L2G_CLOUD_filename)
        L2CLOUDGranule.__init__(self, product_filename=L2G_CLOUD_filename)

    @property
    def tiled_output_variables(self):
        return self.variables

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

    @property
    def QC(self) -> Raster:
        return self.variable(
            variable_name="QC",
            apply_scale=False,
            apply_cloud=False
        ).astype(np.uint16)

    def variable(
            self,
            variable_name: str,
            apply_scale: bool = True,
            apply_cloud: bool = True,
            geometry: RasterGeometry = None,
            **kwargs) -> Raster:

        if variable_name in ["cloud", "water", "QC"]:
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

        if apply_cloud and "cloud" in self.variables:
            cloud = self.variable(
                variable_name="cloud",
                apply_scale=False,
                apply_cloud=False,
                geometry=geometry
            )

            image = rt.where(cloud, np.nan, image)

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

    @classmethod
    def from_swath(
            cls,
            swath_granule: L2CLOUD,
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

        os.makedirs(dirname(output_filename), exist_ok=True)
        dataset_metadata = {}
        output_filename_partial = f"{output_filename}.partial"
        output_filename_rewritten = f"{output_filename}.rewritten"

        if exists(output_filename_partial):
            os.remove(output_filename_partial)

        if variables is None:
            variables = swath_granule.variables

        logger.info(f"processing variables: {', '.join(variables)}")

        with he5py.File(output_filename_partial, "w") as output_file:
            file_timer = Timer()
            logger.info(f"creating HDF-EOS5 file: {cl.file(output_filename_partial)}")

            output_grid = output_file.create_grid(grid_name=grid_name, geometry=gridded_geometry)

            for dataset_name in variables:
                logger.info(f"ingesting dataset: {cl.name(dataset_name)}")
                dataset_timer = Timer()

                source_data = swath_granule.variable(
                    dataset_name,
                    apply_scale=False,
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
                add_offset = 0
                scale_factor = 1

                # if apply_scale:
                #     add_offset = 0
                #     scale_factor = 1
                #
                #     if scaleoffset is not None:
                #         logger.info(f"HDF5 scaleoffset: {scaleoffset}")
                #         kwargs["scaleoffset"] = scaleoffset
                # else:
                #     add_offset = source_metadata["offset"]
                #     scale_factor = source_metadata["scale"]

                if str(source_data.dtype).startswith("float") and scaleoffset is not None:
                    logger.info(f"HDF5 scaleoffset: {scaleoffset}")
                    kwargs["scaleoffset"] = scaleoffset

                logger.info(f"finished ingesting dataset: {cl.name(dataset_name)} ({cl.time(dataset_timer)})")
                write_timer = Timer()
                output_dataset_name = cls.translate_dataset_name(dataset_name)

                if source_data.dtype in (np.float32, np.float64):
                    source_data = source_data.astype(np.float32)
                    _Fillvalue = np.nan

                    logger.info(f"creating dataset: {cl.name(output_dataset_name)} ({cl.val(source_data.dtype)})")

                    output_grid.write_float(
                        field_name=output_dataset_name,
                        image=source_data
                    )
                elif source_data.dtype == np.uint8:
                    logger.info(f"creating dataset: {cl.name(output_dataset_name)} ({cl.val(source_data.dtype)})")

                    if _Fillvalue is None:
                        _Fillvalue = cls.UINT8_FILL

                    output_grid.write_uint8(
                        field_name=output_dataset_name,
                        image=source_data
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
            standard_metadata["AncillaryInputPointer"][...] = "AncillaryNWP"
            orbit = swath_granule.orbit
            scene = swath_granule.scene
            standard_metadata["RegionID"][...] = f"{orbit:05d}_{scene:03d}"
            short_name = L2G_CLOUD_SHORT_NAME
            logger.info(f"L2G CLOUD short name: {cl.val(short_name)}")
            standard_metadata["ShortName"][...] = short_name
            long_name = L2G_CLOUD_LONG_NAME
            logger.info(f"L2G CLOUD long name: {cl.val(long_name)}")
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
            logger.info(f"ImageLineSpacing: {cl.val(cell_size)}")
            standard_metadata["ImageLineSpacing"][...] = cell_size
            logger.info(f"ImagePixelSpacing: {cl.val(cell_size)}")
            standard_metadata["ImagePixelSpacing"][...] = cell_size
            logger.info(f"ImageLines: {cl.val(rows)}")
            standard_metadata["ImageLines"][...] = rows
            logger.info(f"ImagePixels: {cl.val(cols)}")
            standard_metadata["ImagePixels"][...] = cols

            if input_filenames is not None:
                input_pointer = ",".join([basename(filename) for filename in input_filenames])
                logger.info(f"InputPointer: {cl.place(input_pointer)}")
                standard_metadata["InputPointer"][...] = input_pointer

            local_granule_ID = basename(output_filename)
            logger.info(f"LocalGranuleID: {cl.file(local_granule_ID)}")
            standard_metadata["LocalGranuleID"][...] = str(local_granule_ID)
            PROCESSING_LEVEL = "L2G"
            logger.info(f"ProcessingLevelID: {cl.file(PROCESSING_LEVEL)}")
            standard_metadata["ProcessingLevelID"][...] = str(PROCESSING_LEVEL)
            output_group = output_file.data_group(grid_name)

            for dataset_name in variables:
                logger.info(f"ingesting dataset: {cl.name(dataset_name)}")
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
        # shutil.move(output_filename_partial, posterior_filename)
        logger.info(f"completed HDF-EOS5 generation: {cl.file(output_filename)}")
        granule = cls(output_filename)
        XML_metadata_filename = output_filename.replace(".h5", ".h5.met")
        logger.info(f"writing XML metadata file: {cl.file(XML_metadata_filename)}")
        write_XML_metadata(granule.standard_metadata, XML_metadata_filename)

        return granule
