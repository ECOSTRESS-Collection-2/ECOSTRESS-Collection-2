import logging

import numpy as np
from affine import Affine

from MODLAND import find_MODLAND_tiles
import rasters
from rasters import RasterGrid, Raster, RasterGeometry

SINUSOIDAL_PROJECTION = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
MODLAND_TILE_SIZES = {
    250: 4800,
    500: 2400,
    1000: 1200,
    5000: 240,
    10000: 120
}


def calculate_MODLAND_affine(h, v, tile_size):
    # boundaries of sinusodial projection
    UPPER_LEFT_X_METERS = -20015109.355798
    UPPER_LEFT_Y_METERS = 10007554.677899
    LOWER_RIGHT_X_METERS = 20015109.355798
    LOWER_RIGHT_Y_METERS = -10007554.677899

    # size across (width or height) of any equal-area sinusoidal target
    TILE_SIZE_METERS = 1111950.5197665554

    # boundaries of MODIS land grid
    TOTAL_ROWS = 18
    TOTAL_COLUMNS = 36

    y_max = LOWER_RIGHT_Y_METERS + (TOTAL_ROWS - v) * TILE_SIZE_METERS
    x_min = UPPER_LEFT_X_METERS + int(h) * TILE_SIZE_METERS

    cell_size = TILE_SIZE_METERS / tile_size

    # width of pixel
    a = cell_size
    # row rotation
    b = 0.0
    # x-coordinate of upper-left corner of upper-left pixel
    c = x_min
    # column rotation
    d = 0.0
    # height of pixel
    e = -1.0 * cell_size
    # y-coordinate of the upper-left corner of upper-left pixel
    f = y_max
    affine = Affine(a, b, c, d, e, f)

    return affine


def calculate_global_MODLAND_affine(spatial_resolution):
    tile_size = MODLAND_TILE_SIZES[spatial_resolution]
    affine = calculate_MODLAND_affine(0, 0, tile_size)

    return affine


def generate_MODLAND_grid(h, v, tile_size):
    affine = calculate_MODLAND_affine(h, v, tile_size)
    grid = RasterGrid.from_affine(affine, tile_size, tile_size, crs=SINUSOIDAL_PROJECTION)

    return grid


def parsehv(region_name):
    h = int(region_name[1:3])
    v = int(region_name[4:6])

    return h, v


def calculate_global_MODLAND_indices(h, v, tile_size):
    global_modis_affine = calculate_MODLAND_affine(0, 0, tile_size)
    affine = calculate_MODLAND_affine(h, v, tile_size)
    tile_col_indices, tile_row_indices = np.meshgrid(np.arange(tile_size), np.arange(tile_size))
    x, y = affine * (tile_col_indices, tile_row_indices)
    global_col_indices, global_row_indices = ~global_modis_affine * (x, y)

    return global_row_indices, global_col_indices


def calculate_global_MODLAND_serial_indices_hv(h, v, tile_size):
    global_row_indices, global_col_indices = calculate_global_MODLAND_indices(h, v, tile_size)
    serial_index = global_row_indices * tile_size * 36 + global_col_indices
    serial_index = np.int32(serial_index)
    grid = generate_MODLAND_grid(h, v, serial_index.shape[0])
    raster = Raster(serial_index, geometry=grid)

    return raster


def calculate_global_MODLAND_serial_indices(tile, tile_size):
    h = int(tile[1:3])
    v = int(tile[4:6])
    indices = calculate_global_MODLAND_serial_indices_hv(h, v, tile_size)

    return indices


def generate_MODLAND_500m(fine_geometry: RasterGeometry, coarse_resolution: int = None) -> Raster:
    logger = logging.getLogger(__name__)

    FILL_VALUE = -9999
    DEFAULT_COARSE_RESOLUTION = 490

    if coarse_resolution is None:
        coarse_resolution = DEFAULT_COARSE_RESOLUTION

    tiles = find_MODLAND_tiles(fine_geometry.boundary_latlon, return_names=True)
    geometry500m = fine_geometry.rescale(coarse_resolution)
    indices500m = np.full(geometry500m.shape, FILL_VALUE, dtype=np.int32)

    for tile in tiles:
        native_indices = calculate_global_MODLAND_serial_indices(tile, 2400)
        coarse_indices = native_indices.to_geometry(geometry500m, nodata=FILL_VALUE)
        indices500m = rasters.where(indices500m == FILL_VALUE, coarse_indices, indices500m)

    if np.any(indices500m == FILL_VALUE):
        unfilled_proportion = np.count_nonzero(indices500m == FILL_VALUE) / indices500m.size
        logger.warning(f"{round(unfilled_proportion * 100):0.2f} of fine pixels left unfilled")

    indices500m.nodata = FILL_VALUE

    return indices500m


def generate_MODLAND_1000m(fine_geometry: RasterGeometry, coarse_resolution: int = None) -> Raster:
    logger = logging.getLogger(__name__)

    FILL_VALUE = -9999
    DEFAULT_COARSE_RESOLUTION = 980

    if coarse_resolution is None:
        coarse_resolution = DEFAULT_COARSE_RESOLUTION

    tiles = find_MODLAND_tiles(fine_geometry.boundary_latlon, return_names=True)
    geometry = fine_geometry.rescale(coarse_resolution)
    indices = np.full(geometry.shape, FILL_VALUE, dtype=np.int32)

    for tile in tiles:
        native_indices = calculate_global_MODLAND_serial_indices(tile, 1200)
        coarse_indices = native_indices.to_geometry(geometry, nodata=FILL_VALUE)
        indices = rasters.where(indices == FILL_VALUE, coarse_indices, indices)

    if np.any(indices == FILL_VALUE):
        unfilled_proportion = np.count_nonzero(indices == FILL_VALUE) / indices.size
        logger.warning(f"{round(unfilled_proportion * 100):0.2f} of fine pixels left unfilled")

    indices.nodata = FILL_VALUE

    return indices


def calculate_global_MODLAND_columns(spatial_resolution):
    tile_size = MODLAND_TILE_SIZES[spatial_resolution]
    global_MODLAND_columns = tile_size * 36

    return global_MODLAND_columns


def generate_MODLAND_indices(geometry: RasterGeometry, spatial_resolution: float) -> Raster:
    x, y = geometry.get_xy(projection=SINUSOIDAL_PROJECTION)
    global_MODLAND_affine = calculate_global_MODLAND_affine(spatial_resolution)
    global_col_indices, global_row_indices = ~global_MODLAND_affine * (x, y)
    global_col_indices = global_col_indices.astype(np.int32)
    global_row_indices = global_row_indices.astype(np.int32)
    global_MODLAND_columns = calculate_global_MODLAND_columns(spatial_resolution)
    serial_index = global_row_indices * global_MODLAND_columns + global_col_indices
    serial_index = serial_index.astype(np.int32)
    serial_index = Raster(serial_index, geometry=geometry)

    return serial_index
