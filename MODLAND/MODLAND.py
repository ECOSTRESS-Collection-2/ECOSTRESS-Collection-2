"""
This module handles georeferencing MODIS land tiles.
"""
# rotation in affine transformation is not supported
# latitude and longitude are in WGS84 geographic coordinate system
import warnings
import geopandas as gpd
import numpy as np
import pyproj
import rasterio
from affine import Affine
from numpy import where, nan, isnan, all
from pyproj import Proj
from scipy.spatial.qhull import ConvexHull
from shapely.geometry import Polygon, Point, MultiPolygon
from shapely.geometry import mapping
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import LinearRing

# from transform import transform_shape

__author__ = 'Gregory H. Halverson'

CENTER_PIXEL_COORDINATES_DEFAULT = True

MODIS_LAND_TILE_PROJECTION_WKT = 'PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not_specified_based_on_custom_spheroid",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
MODIS_LAND_TILE_PROJECTION_PROJ4 = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
MODIS_LAND_TILE_PROJECTION_PCI = ['SIN         E700', 'METRE',
                                  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]
MODIS_LAND_TILE_PROJECTION_MI = 'Earth Projection 16, 104, "m", 0'
MODIS_LAND_TILE_PROJECTION_EPSG = 6842

SINUSOIDAL = Proj(MODIS_LAND_TILE_PROJECTION_PROJ4)
WGS84 = Proj('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ')


def get_proj4(projection: Proj or str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if isinstance(projection, Proj):
            result = projection.crs.to_proj4()

        elif isinstance(projection, rasterio.crs.CRS):
            result = projection.to_proj4()

        elif isinstance(projection, pyproj.crs.crs.CRS):
            result = projection.to_proj4()

        elif isinstance(projection, str):
            # result = pycrs.parse.from_unknown_text(projection).to_proj4()
            result = pyproj.CRS(projection).to_proj4()

        else:
            raise ValueError(f"projection not recognized ({type(projection)})")

    result = result.replace("+shell=epsg:", "epsg:")

    return result


def transform_shape(
        shape: BaseGeometry,
        source_projection: Proj or str,
        target_projection: Proj or str) -> BaseGeometry:
    if not isinstance(shape, BaseGeometry):
        raise ValueError("invalid shape")

    # TODO need to stop relying on deprecated proj4 as common projection encoding
    source_projection = get_proj4(source_projection)
    target_projection = get_proj4(target_projection)
    projected_shape = gpd.GeoDataFrame({}, geometry=[shape], crs=source_projection).to_crs(target_projection).geometry[
        0]

    return projected_shape


def transform_sinusoidal_to_latlon(shape):
    return transform_shape(shape, SINUSOIDAL, WGS84)


def transform_latlon_to_sinusoidal(shape):
    return transform_shape(shape, WGS84, SINUSOIDAL)


# constants and formulas based on
# https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
# with slight modifications to boundaries based on proj4

# spherical earth radius (authalic radius)
EARTH_RADIUS_METERS = 6371007.181000

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


# transforms WGS84 latitude and longitude to sinusoidal coordinates in meters
def latlon_to_sinusoidal(latitude, longitude):
    if latitude < -90 or latitude > 90:
        raise ValueError('latitude (%f) out of bounds' % latitude)

    if longitude < -180 or longitude > 180:
        raise ValueError('longitude (%f) out of bounds' % longitude)

    return SINUSOIDAL(longitude, latitude)


# transforms sinusoidal coordinates in meters to WGS84 latitude and longitude
def sinusoidal_to_latlon(x_sinusoidal, y_sinusoidal):
    if x_sinusoidal < UPPER_LEFT_X_METERS or x_sinusoidal > LOWER_RIGHT_X_METERS:
        raise ValueError('sinusoidal x coordinate (%f) out of bounds' % x_sinusoidal)

    if y_sinusoidal < LOWER_RIGHT_Y_METERS or y_sinusoidal > UPPER_LEFT_Y_METERS:
        raise ValueError('sinusoidal y (%f) coordinate out of bounds' % y_sinusoidal)

    longitude, latitude = SINUSOIDAL(x_sinusoidal, y_sinusoidal, inverse=True)

    if x_sinusoidal < 0 and x_sinusoidal < SINUSOIDAL(-180, latitude)[0]:
        return None, None

    if x_sinusoidal > 0 and x_sinusoidal > SINUSOIDAL(180, latitude)[0]:
        return None, None

    return latitude, longitude


# MODIS land target indices for target containing sinusoidal coordinate
def sinusoidal_to_MODLAND(x_sinusoidal, y_sinusoidal):
    if x_sinusoidal < UPPER_LEFT_X_METERS or x_sinusoidal > LOWER_RIGHT_X_METERS:
        raise ValueError('sinusoidal x coordinate (%f) out of bounds' % x_sinusoidal)

    if y_sinusoidal < LOWER_RIGHT_Y_METERS or y_sinusoidal > UPPER_LEFT_Y_METERS:
        raise ValueError('sinusoidal y (%f) coordinate out of bounds' % y_sinusoidal)

    horizontal_index = int((x_sinusoidal - UPPER_LEFT_X_METERS) / TILE_SIZE_METERS)
    vertical_index = int((-1 * (y_sinusoidal + LOWER_RIGHT_Y_METERS)) / TILE_SIZE_METERS)

    if horizontal_index == TOTAL_COLUMNS:
        horizontal_index -= 1

    if vertical_index == TOTAL_ROWS:
        vertical_index -= 1

    return horizontal_index, vertical_index


# MODIS land target indices for target containing latitude and longitude
def latlon_to_MODLAND(latitude, longitude):
    if latitude < -90 or latitude > 90:
        raise ValueError('latitude (%f) out of bounds' % latitude)

    if longitude < -180 or longitude > 180:
        raise ValueError('longitude (%f) out of bounds' % longitude)

    return sinusoidal_to_MODLAND(*latlon_to_sinusoidal(latitude, longitude))


# x coordinate of left side of sinusoidal target
# (upper-left corner of pixel)
def MODLAND_left_x_meters(horizontal_index):
    if horizontal_index >= TOTAL_COLUMNS or horizontal_index < 0:
        raise IndexError('horizontal index (%d) out of bounds' % horizontal_index)

    return UPPER_LEFT_X_METERS + int(horizontal_index) * TILE_SIZE_METERS


# x coordinate of right side of right-most pixels of sinusoidal target
def MODLAND_right_x_meters(horizontal_index):
    if horizontal_index >= TOTAL_COLUMNS or horizontal_index < 0:
        raise IndexError('horizontal (%d) index out of bounds' % horizontal_index)

    return MODLAND_left_x_meters(horizontal_index) + TILE_SIZE_METERS


# y coordinate of top side of sinusoidal target
# (upper-left corner of pixel)
def MODLAND_top_y_meters(vertical_index):
    if vertical_index >= TOTAL_ROWS or vertical_index < 0:
        raise IndexError('vertical index (%d) out of bounds' % vertical_index)

    return LOWER_RIGHT_Y_METERS + (TOTAL_ROWS - vertical_index) * TILE_SIZE_METERS


# y coordinate of the bottom side of bottom-most pixels of sinusoidal target
def MODLAND_bottom_y_meters(vertical_index):
    if vertical_index >= TOTAL_ROWS or vertical_index < 0:
        raise IndexError('vertical index (%d) out of bounds' % vertical_index)

    return LOWER_RIGHT_Y_METERS + (TOTAL_ROWS - 1 - vertical_index) * TILE_SIZE_METERS


# size across each cell in meters given the number of cells across the target
def MODLAND_cell_size_meters(cells_across_tile):
    return TILE_SIZE_METERS / cells_across_tile


# encapsulation of MODIS land target at given indices
# affine transform of raster can be calculated given count of rows and columns
class MODISLandTile:
    def __init__(self, horizontal_index, vertical_index, rows, columns):

        if horizontal_index >= TOTAL_COLUMNS or horizontal_index < 0:
            raise IndexError('horizontal index (%d) out of bounds' % horizontal_index)

        if vertical_index >= TOTAL_ROWS or vertical_index < 0:
            raise IndexError('vertical index (%d) out of bounds' % vertical_index)

        if rows < 0:
            raise ValueError('rows cannot be negative (%d)' % rows)

        if columns < 0:
            raise ValueError('columns cannot be negative (%d)' % columns)

        self.horizontal_index = horizontal_index
        self.vertical_index = vertical_index
        self.rows = rows
        self.columns = columns

    def __str__(self):
        return "<MODISLandTile h%02dv%02d, rows: %d, columns: %d>" \
               % (self.horizontal_index, self.vertical_index, self.rows, self.columns)

    # x coordinate of left side of sinusoidal target
    # (upper-left corner of pixel)
    @property
    def x_min(self):
        return MODLAND_left_x_meters(self.horizontal_index)

    # x coordinate of center of upper-left pixel
    @property
    def x_min_center(self):
        return self.x_min + self.cell_width_meters / 2.0

    # x coordinate of right side right-most pixels of sinusoidal target
    @property
    def x_max(self):
        return MODLAND_right_x_meters(self.horizontal_index)

    # y coordinate of top side of sinusoidal target
    # (upper-left corner of pixel)
    @property
    def y_max(self):
        return MODLAND_top_y_meters(self.vertical_index)

    # y coordinate of center of upper-left pixel
    @property
    def y_max_center(self):
        return self.y_max - self.cell_height_meters / 2.0

    # y coordinate of the bottom side of bottom-most pixels of sinusoidal target
    @property
    def y_min(self):
        return MODLAND_bottom_y_meters(self.vertical_index)

    # width of cell in meters given number of columns
    @property
    def cell_width_meters(self):
        return MODLAND_cell_size_meters(self.columns)

    # positive height of cell in meters given number of rows
    @property
    def cell_height_meters(self):
        return MODLAND_cell_size_meters(self.rows)

    # tuple of cell width and height
    @property
    def cell_size_meters(self):
        return (self.cell_width_meters, self.cell_height_meters)

    # affine transform (as tuple) of target given cells across target
    @property
    def affine_tuple(self):

        # width of pixel
        a = self.cell_width_meters

        # row rotation
        b = 0.0

        # x-coordinate of upper-left corner of upper-left pixel
        c = self.x_min

        # column rotation
        d = 0.0

        # height of pixel
        e = -1.0 * self.cell_height_meters

        # y-coordinate of the upper-left corner of upper-left pixel
        f = self.y_max

        affine_transform = (a, b, c, d, e, f)

        return affine_transform

    # affine transform as affine.Affine
    @property
    def affine_transform(self):
        return Affine(*self.affine_tuple)

    @property
    def affine_center(self):
        return self.affine_transform * Affine.translation(0.5, 0.5)

    # affine transform as tuple in ESRI world format
    @property
    def affine_esri(self):

        # width of pixel
        a = self.cell_width_meters

        # column rotation
        d = 0.0

        # row rotation
        b = 0.0

        # height of pixel
        e = -1.0 * self.cell_height_meters

        # x coordinate of center of upper-left pixel
        c = self.x_min_center

        # y coordinate of center of upper-left pixel
        f = self.y_max_center

        affine_transform = (a, d, b, e, c, f)

        return affine_transform

    # save esri world file for target
    def save_world_file(self, filename):
        with open(filename, 'w') as f:
            f.write('\n'.join([str(parameter) for parameter in (self.affine_esri)]))

    # affine transform as tuple in GDAL format
    @property
    def affine_gdal(self):

        # x-coordinate of upper-left corner of upper-left pixel
        c = self.x_min

        # width of pixel
        a = self.cell_width_meters

        # row rotation
        b = 0.0

        # y-coordinate of the upper-left corner of upper-left pixel
        f = self.y_max

        # column rotation
        d = 0.0

        # height of pixel
        e = -1.0 * self.cell_height_meters

        affine_transform = (c, a, b, f, d, e)

        return affine_transform

    # x, y sinusoidal coordinate in meters from row and column
    # top left of pixel if center is set to False
    # center of pixel if center is set to True
    def sinusoidal(self, row, column, center=CENTER_PIXEL_COORDINATES_DEFAULT):
        if row >= self.rows or row < 0:
            raise IndexError('row (%d) out of bounds' % row)

        if column >= self.columns or column < 0:
            raise IndexError('column (%d) out of bounds' % column)

        # top-left corner coordinates of pixel
        x_sinusoidal, y_sinusoidal = self.affine_transform * (column, row)

        # offset coordinates to center of pixel
        if center:
            x_sinusoidal += self.cell_width_meters / 2.0
            y_sinusoidal -= self.cell_height_meters / 2.0

        return (x_sinusoidal, y_sinusoidal)

    # get row and column of cell at sinusoidal coordinates
    def row_column_from_sinusoidal(self, x_sinusoidal, y_sinusoidal):
        if x_sinusoidal < self.x_min or x_sinusoidal >= self.x_max:
            raise ValueError('sinusoidal x coordinate (%f) out of bounds' % x_sinusoidal)

        if y_sinusoidal <= self.y_min or y_sinusoidal > self.y_max:
            raise ValueError('sinusoidal y (%f) coordinate out of bounds' % y_sinusoidal)

        return tuple([int(index)
                      for index
                      in ~(self.affine_transform) * (x_sinusoidal, y_sinusoidal)])

    def row_column_from_latlon(self, latitude, longitude):
        return self.row_column_from_sinusoidal(*latlon_to_sinusoidal(latitude, longitude))

    # get latitude and longitude of cell at row and column
    # top left of pixel if center is set to False
    # center of pixel if center is set to True
    def latlon(self, row, column, center=CENTER_PIXEL_COORDINATES_DEFAULT):
        if row >= self.rows or row < 0:
            raise IndexError('row (%d) out of bounds' % row)

        if column >= self.columns or column < 0:
            raise IndexError('column (%d) out of bounds' % column)

        latitude, longitude = sinusoidal_to_latlon(*self.sinusoidal(row, column, center=center))

        return (latitude, longitude)

    # latitude matrix and longitude matrix for coordinates of each cell in target
    # top left of pixel if center is set to False
    # center of pixel if center is set to True
    def latlon_matrices(self, center=CENTER_PIXEL_COORDINATES_DEFAULT):
        # lon, lat = sinusoidal_projection(*self.sinusoidal_matrices(center=center), inverse=True)
        _, _, lon, lat = self.sinusoidal_matrices(center=center)

        return lat, lon

    # sinusoidal x and y matrices of each cell in target
    # top left of pixel if center is set to False
    # center of pixel if center is set to True
    def sinusoidal_matrices(self, center=CENTER_PIXEL_COORDINATES_DEFAULT):

        if center:
            affine = self.affine_center
        else:
            affine = self.affine_transform

        x_matrix, y_matrix = np.meshgrid(np.arange(self.columns), np.arange(self.rows)) * affine
        lon, lat = SINUSOIDAL(x_matrix, y_matrix, inverse=True)

        if self.horizontal_index < 18:
            valid = lon < 0
        else:
            valid = lon >= 0

        x_matrix = where(valid, x_matrix, nan)
        y_matrix = where(valid, y_matrix, nan)
        lon = where(valid, lon, nan)
        lat = where(valid, lat, nan)

        return x_matrix, y_matrix, lon, lat

    # checks if a sinusoidal coordinate falls within the target
    def contains_sinusoidal(self, x_sinusoidal, y_sinusoidal):
        if x_sinusoidal < UPPER_LEFT_X_METERS or x_sinusoidal > LOWER_RIGHT_X_METERS:
            raise ValueError('sinusoidal x coordinate (%f) out of bounds' % x_sinusoidal)

        if y_sinusoidal < LOWER_RIGHT_Y_METERS or y_sinusoidal > UPPER_LEFT_Y_METERS:
            raise ValueError('sinusoidal y (%f) coordinate out of bounds' % y_sinusoidal)

        h, v = sinusoidal_to_MODLAND(x_sinusoidal, y_sinusoidal)

        return h == self.horizontal_index and v == self.vertical_index

    # checks if a latitude and longitude coordinate falls within the target
    def contains_latlon(self, latitude, longitude):
        return self.contains_sinusoidal(*latlon_to_sinusoidal(latitude, longitude))

    @property
    def outline_sinusoidal(self):
        upper_left_x = lower_left_x = self.x_min
        upper_right_x = lower_right_x = self.x_max
        upper_left_y = upper_right_y = self.y_max
        lower_left_y = lower_right_y = self.y_min

        # tests if corners exist
        upper_left_valid = any(sinusoidal_to_latlon(upper_left_x, upper_left_y))
        upper_right_valid = any(sinusoidal_to_latlon(upper_right_x, upper_right_y))
        lower_right_valid = any(sinusoidal_to_latlon(lower_right_x, lower_right_y))
        lower_left_valid = any(sinusoidal_to_latlon(lower_left_x, lower_left_y))

        # case where all corners exist
        if all([upper_left_valid, upper_right_valid, lower_right_valid, lower_left_valid]):
            return Polygon(LinearRing([
                [upper_left_x, upper_left_y],
                [upper_right_x, upper_right_y],
                [lower_right_x, lower_right_y],
                [lower_left_x, lower_left_y]
            ]))
        else:
            x_matrix, y_matrix, lon, lat = self.sinusoidal_matrices()
            points = np.dstack([x_matrix.flatten(), y_matrix.flatten()])[0]
            points = points[~all(isnan(points), axis=1)]
            hull = ConvexHull(points)
            vertices = hull.points[hull.vertices]
            poly = Polygon(vertices)
            return poly

    @property
    def outline_latlon(self):
        return transform_sinusoidal_to_latlon(self.outline_sinusoidal)

    @property
    def bounds_sinusoidal(self):
        return self.x_min, self.y_min, self.x_max, self.y_max


# finds target containing a latitude and longitude coordinate
# then finds pixel within that target nearest to coordinate given size of matrix
# returns horizontal index, vertical index, row, and column as tuple
def latlon_to_MODLAND_pixel(latitude, longitude, rows_per_tile, columns_per_tile):
    horizontal_index, vertical_index = latlon_to_MODLAND(latitude, longitude)
    tile = MODISLandTile(horizontal_index, vertical_index, rows_per_tile, columns_per_tile)
    row, column = tile.row_column_from_latlon(latitude, longitude)

    return horizontal_index, vertical_index, row, column


# finds target containing a sinusoidal coordinate
# then finds pixel within that target nearest to coordinate given size of matrix
# returns horizontal index, vertical index, row, and column as tuple
def sinusoidal_to_MODLAND_pixel(x_sinusoidal, y_sinusoidal, rows_per_tile, columns_per_tile):
    horizontal_index, vertical_index = sinusoidal_to_MODLAND(x_sinusoidal, y_sinusoidal)
    tile = MODISLandTile(horizontal_index, vertical_index, rows_per_tile, columns_per_tile)
    row, column = tile.row_column_from_sinusoidal(x_sinusoidal, y_sinusoidal)

    return horizontal_index, vertical_index, row, column


# calculate polygon outline of MODIS land target in sinusoidal coordinates
def outline_sinusoidal(h, v):
    """
    Calculate polygon outline of MODIS land target in sinusoidal coordinates.
    :param h: horizontal index of MODIS land target
    :param v: vertical index of MODIS land target
    :return: shapely geometry object containing polygon outline of MODIS land target in sinusoidal coordinates
    """
    return MODISLandTile(h, v, 100, 100).outline_sinusoidal()


# calculate polygon outline of MODIS land target in latitude and longitude
def outline_latlon(h, v):
    """
    Calculate polygon outline of MODIS land target in latitude and longitude.
    :param h: horizontal index of MODIS land target
    :param v: vertical index of MODIS land target
    :return: shapely geometry object containing polygon outline of MODIS land target in latitude and longitude
    """
    return MODISLandTile(h, v, 100, 100).outline_latlon


# calculate MODIS land tiles intersecting a polygon in latitude and longitude
def find_MODLAND_tiles(target_geometry_latlon, return_names=True, land_only=True):
    """
    Calculate MODIS land tiles intersecting a polygon in latitude and longitude.
    :param target_geometry_latlon: target polygon as a shapely geometry object with latitude and longitude coordinates
    :return: set of tuples of h and v indices of MODIS land tiles intersecting target polygon
    """

    # list of tiles intersecting target polygon
    intersecting_tiles = []

    # list of tiles at boundary coordinates of target polygon
    boundary_tiles = []

    # calculate target at each point in the target polygon

    if isinstance(target_geometry_latlon, Polygon):
        shapes = mapping(target_geometry_latlon)['coordinates']
    elif isinstance(target_geometry_latlon, MultiPolygon):
        shapes = mapping(target_geometry_latlon.convex_hull)['coordinates']
    elif isinstance(target_geometry_latlon, Point):
        shapes = [target_geometry_latlon]
    else:
        raise ValueError("unsupported target geometry: {}".format(target_geometry_latlon))

    # iterate through shapes
    for shape in shapes:

        if isinstance(shape, Point):
            shape = [shape]

        # iterate through coordinates
        for coordinate in shape:
            if isinstance(coordinate, Point):
                longitude = coordinate.x
                latitude = coordinate.y
            else:
                # pull latitude and longitude from coordinate
                longitude, latitude = coordinate

            # calculate target at coordinate
            tile = latlon_to_MODLAND(latitude, longitude)

            # add target to list
            boundary_tiles += [tile]

    # set of tiles at boundary coordinates of target polygon
    boundary_tiles = set(boundary_tiles)

    # pull h and v indices from set of tiles
    horizontal_indices, vertical_indices = zip(*boundary_tiles)

    # iterate through range of horizontal indices
    for h in range(min(horizontal_indices), max(horizontal_indices) + 1):

        # iterate through range of vertical indices
        for v in range(min(vertical_indices), max(vertical_indices) + 1):
            polygon = outline_latlon(h, v)

            # check if target intersects target polygon
            if polygon.intersects(target_geometry_latlon):
                # add target to list
                intersecting_tiles += [(h, v)]

    # set of tiles intersecting target polygon
    intersecting_tiles = set(intersecting_tiles)

    if return_names:
        return sorted(["h{:02d}v{:02d}".format(h, v) for h, v in intersecting_tiles])
    else:
        return intersecting_tiles
