from typing import List, Any, Union
from os.path import exists, join
from glob import glob
import numpy as np
import rasters as rt
from rasters import Raster, RasterGrid, RasterGeolocation, KDTree
from timer import Timer
import logging
import colored_logging as cl

logger = logging.getLogger(__name__)


def calculate_scan_start(scan_number: int) -> int:
    """
    Calculate zero-based row index for zero-based scan index
    """
    return scan_number * 128

def calculate_scan_stop(scan_number: int) -> int:
    """
    Calculate zero-based exclusive limit row index for zero-based scan index
    """
    return 128 * (scan_number + 1)

def subset_scan_from_swath(swath_raster: Raster, scan_number: int) -> Raster:
    """
    Subset swath data layer to zero-based scan index
    """
    return swath_raster[calculate_scan_start(scan_number):calculate_scan_stop(scan_number), :]

def subset_head_from_scan(scan_raster: Raster) -> Raster:
    """
    Subset 28-row overlapping head from scan
    """
    return scan_raster[:28, :]

def subset_nadir_from_scan(scan_raster: Raster) -> Raster:
    """
    Subset 72-row nadir from scan
    """
    return scan_raster[28:100, :]

def subset_tail_from_scan(scan_raster: Raster) -> Raster:
    """
    Subset 28-row overlapping tail from scan
    """
    return scan_raster[100:, :]


def generate_scan_kd_trees(swath_geometry: RasterGeolocation, cell_size_degrees: float = 0.0007) -> List[KDTree]:
    scan_trees = []

    for scan in range(44):
        logger.info(f"started generating K-D tree for scan {scan} ({scan + 1} / 44)")
        timer = Timer()
        start_index = 128 * scan
        end_index = 128 * (scan + 1)
        scan_swath_geometry = swath_geometry[start_index:end_index, :]
        scan_gridded_geometry = scan_swath_geometry.geographic(cell_size_degrees)
        scan_tree = KDTree(scan_swath_geometry, scan_gridded_geometry)
        scan_trees.append(scan_tree)
        logger.info(f"finished generating K-D tree for scan {scan} ({scan + 1} / 44) ({timer})")

    return scan_trees

def clip_tails(swath: Union[Raster, RasterGeolocation, np.ndarray]):
    if isinstance(swath, Raster):
        return swath.contain(clip_tails(swath.array), geometry=clip_tails(swath.geometry))
    elif isinstance(swath, RasterGeolocation):
        return RasterGeolocation(x=clip_tails(swath.lon), y=clip_tails(swath.lat))
    elif isinstance(swath, np.ndarray):
        return swath[np.arange(5632) % 128 < 105, :]
    else:
        raise ValueError(f"invalid swath {type(swath)}")

def resample_scan_by_scan(
        swath_raster: Raster,
        target_geometry: RasterGrid = None,
        scan_kd_trees: Union[List[KDTree], str] = None,  # either a list of K-D tree object or a directory name
        merge_method: str = None,
        nodata: Any = np.nan,
        cell_size_degrees: float = 0.0007,
        search_radius_meters: float = 100) -> Raster:
    gridded_image = None
    swath_geometry = swath_raster.geometry
    dtype = swath_raster.dtype

    if nodata is np.nan and "int" in str(swath_raster.dtype):
        raise ValueError("cannot use NaN as nodata value for integer layer")

    if target_geometry is None:
        target_geometry = swath_geometry.geographic(cell_size_degrees)

    if merge_method is None:
        if "int" in str(dtype):
            merge_method = "or"
            gridded_image = Raster(np.full(target_geometry.shape, 0, dtype=dtype), geometry=target_geometry)

            if nodata is np.nan:
                if "uint8" in str(dtype):
                    nodata = 255
                elif "uint16" in str(dtype):
                    nodata = 65535
        else:
            merge_method = "average"

    if scan_kd_trees is None:
        logger.warning("scan K-D trees not supplied")
        scan_kd_trees = generate_scan_kd_trees(swath_geometry, cell_size_degrees=cell_size_degrees)
    elif isinstance(scan_kd_trees, str) and exists(scan_kd_trees):
        kd_tree_directory = scan_kd_trees

        scan_kd_trees = [KDTree.load(filename) for filename in sorted(glob(join(kd_tree_directory, "*.kdtree")))]

    for scan in range(44):
        logger.info(f"started resampling scan {scan} ({scan + 1} / 44)")
        timer = Timer()
        start_index = 128 * scan
        end_index = 128 * (scan + 1)
        scan_swath_subset = swath_raster[start_index:end_index, :]
        scan_swath_geometry = swath_geometry[start_index:end_index, :]
        scan_gridded_geometry = scan_swath_geometry.geographic(cell_size_degrees)
        scan_tree = scan_kd_trees[scan]

        scan_gridded_local = scan_swath_subset.to_geometry(
            scan_gridded_geometry,
            kd_tree=scan_tree,
            search_radius_meters=search_radius_meters,
            nodata=nodata
        )

        scan_gridded = scan_gridded_local.to_geometry(target_geometry, search_radius_meters=search_radius_meters)

        if gridded_image is None:
            gridded_image = scan_gridded
        else:
            if merge_method == "or":
                gridded_image = gridded_image | scan_gridded
            else:
                gridded_image = rt.where(
                    np.isnan(gridded_image), scan_gridded, gridded_image)

                overlap_average = (gridded_image + scan_gridded) / 2

                gridded_image = rt.where(
                    np.isnan(overlap_average), gridded_image, overlap_average)

        logger.info(f"finished resampling scan {scan} ({scan + 1} / 44) ({timer})")

    return gridded_image