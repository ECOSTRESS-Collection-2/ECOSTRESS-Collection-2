import base64
import io
import json
import logging
import sys
import warnings
from datetime import datetime
from datetime import timedelta, date
from glob import glob
from math import floor
from os import makedirs
from os.path import abspath, dirname, expanduser
from os.path import basename
from os.path import join, exists
from os.path import splitext
from time import perf_counter
from typing import List, Set, Union, Tuple
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.crs
import shapely
import shapely.geometry
import shapely.wkt
import xmltodict
from affine import Affine
from dateutil import parser
import mgrs
from rasterio.features import rasterize
from rasterio.warp import Resampling
from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt
from rasters import Point, Polygon, BBox
from shapely.geometry.base import BaseGeometry
from six import string_types

import colored_logging as cl
from rasters import RasterGrid, Raster, RasterGeometry, WGS84, CRS

# from transform.UTM import UTM_proj4_from_latlon, UTM_proj4_from_zone

pd.options.mode.chained_assignment = None  # default='warn'

DEFAULT_ALBEDO_RESOLUTION = 10
DEFAULT_SEARCH_DAYS = 10
DEFAULT_CLOUD_MIN = 0
DEFAULT_CLOUD_MAX = 50
DEFAULT_ORDER_BY = "-beginposition"

SENTINEL_POLYGONS_FILENAME = join(abspath(dirname(__file__)), "sentinel2_tiles_world_with_land.geojson")

DEFAULT_WORKING_DIRECTORY = "."
DEFAULT_DOWNLOAD_DIRECTORY = "sentinel_download"
DEFAULT_PRODUCTS_DIRECTORY = "sentinel_products"

logger = logging.getLogger(__name__)

def UTM_proj4_from_latlon(lat: float, lon: float) -> str:
    UTM_zone = (floor((lon + 180) / 6) % 60) + 1
    UTM_proj4 = f"+proj=utm +zone={UTM_zone} {'+south ' if lat < 0 else ''}+ellps=WGS84 +datum=WGS84 +units=m +no_defs"

    return UTM_proj4


def UTM_proj4_from_zone(zone: str):
    zone_number = int(zone[:-1])

    if zone[-1].upper() == "N":
        hemisphere = ""
    elif zone[-1].upper() == "S":
        hemisphere = "+south "
    else:
        raise ValueError(f"invalid hemisphere in zone: {zone}")

    UTM_proj4 = f"+proj=utm +zone={zone_number} {hemisphere}+datum=WGS84 +units=m +no_defs"

    return UTM_proj4


def load_geojson_as_wkt(geojson_filename: str) -> str:
    return geojson_to_wkt(read_geojson(geojson_filename))


def parse_sentinel_granule_id(granule_id: str) -> dict:
    # Compact Naming Convention
    #
    # The compact naming convention is arranged as follows:
    #
    # MMM_MSIXXX_YYYYMMDDHHMMSS_Nxxyy_ROOO_Txxxxx_<Product Discriminator>.SAFE
    #
    # The products contain two dates.
    #
    # The first date (YYYYMMDDHHMMSS) is the datatake sensing time.
    # The second date is the "<Product Discriminator>" field, which is 15 characters in length, and is used to distinguish between different end user products from the same datatake. Depending on the instance, the time in this field can be earlier or slightly later than the datatake sensing time.
    #
    # The other components of the filename are:
    #
    # MMM: is the mission ID(S2A/S2B)
    # MSIXXX: MSIL1C denotes the Level-1C product level/ MSIL2A denotes the Level-2A product level
    # YYYYMMDDHHMMSS: the datatake sensing start time
    # Nxxyy: the PDGS Processing Baseline number (e.g. N0204)
    # ROOO: Relative Orbit number (R001 - R143)
    # Txxxxx: Tile Number field
    # SAFE: Product Format (Standard Archive Format for Europe)
    #
    # Thus, the following filename
    #
    # S2A_MSIL1C_20170105T013442_N0204_R031_T53NMJ_20170105T013443.SAFE
    #
    # Identifies a Level-1C product acquired by Sentinel-2A on the 5th of January, 2017 at 1:34:42 AM. It was acquired over Tile 53NMJ(2) during Relative Orbit 031, and processed with PDGS Processing Baseline 02.04.
    #
    # In addition to the above changes, a a TCI (True Colour Image) in JPEG2000 format is included within the Tile folder of Level-1C products in this format. For more information on the TCI, see the Definitions page here.
    # https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/naming-convention
    parts = granule_id.split("_")

    return {
        "mission_id": parts[0],
        "product": parts[1],
        "date": parser.parse(parts[2]),
        "baseline": parts[3],
        "orbit": parts[4],
        "tile": parts[5][1:]
    }


def resize_affine(affine: Affine, cell_size: float) -> Affine:
    if not isinstance(affine, Affine):
        raise ValueError("invalid affine transform")

    new_affine = Affine(cell_size, affine.b, affine.c, affine.d, -cell_size, affine.f)

    return new_affine


def UTC_to_solar(time_UTC, lon):
    return time_UTC + timedelta(hours=(np.radians(lon) / np.pi * 12))

class MGRS(mgrs.MGRS):
    def bbox(self, tile: str) -> BBox:
        if len(tile) == 5:
            precision = 100000
        elif len(tile) == 7:
            precision = 10000
        elif len(tile) == 9:
            precision = 1000
        elif len(tile) == 11:
            precision = 100
        elif len(tile) == 13:
            precision = 10
        elif len(tile) == 15:
            precision = 1
        else:
            raise ValueError(f"unrecognized MGRS tile: {tile}")

        zone, hemisphere, xmin, ymin = self.MGRSToUTM(tile)
        crs = CRS(UTM_proj4_from_zone(f"{int(zone)}{str(hemisphere).upper()}"))
        xmax = xmin + precision
        ymax = ymin + precision

        bbox = BBox(
            xmin=xmin,
            ymin=ymin,
            xmax=xmax,
            ymax=ymax,
            crs=crs
        )

        return bbox


class SentinelTileGrid(MGRS):
    def __init__(self, *args, target_resolution: float = 30, **kwargs):
        super(SentinelTileGrid, self).__init__()
        self.target_resolution = target_resolution
        self._sentinel_polygons = None

    def __repr__(self) -> str:
        return f"SentinelTileGrid(target_resolution={self.target_resolution})"

    @property
    def sentinel_polygons(self) -> gpd.GeoDataFrame:
        if self._sentinel_polygons is None:
            self._sentinel_polygons = gpd.read_file(SENTINEL_POLYGONS_FILENAME)

        return self._sentinel_polygons

    @property
    def crs(self) -> CRS:
        return CRS(self._sentinel_polygons.crs)

    def UTM_proj4(self, tile: str) -> str:
        zone, hemisphere, _, _ = self.MGRSToUTM(tile)
        proj4 = UTM_proj4_from_zone(f"{int(zone)}{str(hemisphere).upper()}")

        return proj4

    def footprint(
            self,
            tile: str,
            in_UTM: bool = False,
            round_UTM: bool = True,
            in_2d: bool = True) -> Polygon:
        try:
            polygon = Polygon(self.sentinel_polygons[self.sentinel_polygons.Name == tile].iloc[0]["geometry"], crs=self.crs)
        except Exception as e:
            raise ValueError(f"polygon for target {tile} not found")

        if in_2d:
            polygon = Polygon([xy[0:2] for xy in polygon.exterior.coords], crs=self.crs)

        if in_UTM:
            UTM_proj4 = self.UTM_proj4(tile)
            # print(f"transforming: {WGS84} -> {UTM_proj4}")
            polygon = polygon.to_crs(UTM_proj4)

            if round_UTM:
                polygon = Polygon([[round(item) for item in xy] for xy in polygon.exterior.coords], crs=polygon.crs)

        return polygon

    def footprint_UTM(self, tile: str) -> Polygon:
        return self.footprint(
            tile=tile,
            in_UTM=True,
            round_UTM=True,
            in_2d=True
        )

    def bbox(self, tile: str, MGRS: bool = False) -> BBox:
        if len(tile) != 5 or MGRS:
            return super(SentinelTileGrid, self).bbox(tile=tile)

        polygon = self.footprint(
            tile=tile,
            in_UTM=True,
            round_UTM=True,
            in_2d=True
        )

        bbox = polygon.bbox

        return bbox

    def tiles(self, target_geometry: shapely.geometry.shape) -> Set[str]:
        if isinstance(target_geometry, str):
            target_geometry = shapely.wkt.loads(target_geometry)

        matches = self.sentinel_polygons[self.sentinel_polygons.intersects(target_geometry)]
        tiles = set(sorted(list(matches.apply(lambda row: row.Name, axis=1))))

        return tiles

    def tile_footprints(
            self,
            target_geometry: shapely.geometry.shape or gpd.GeoDataFrame,
            calculate_area: bool = False,
            eliminate_redundancy: bool = False) -> gpd.GeoDataFrame:
        if isinstance(target_geometry, str):
            target_geometry = shapely.wkt.loads(target_geometry)

        if isinstance(target_geometry, BaseGeometry):
            target_geometry = gpd.GeoDataFrame(geometry=[target_geometry], crs="EPSG:4326")

        if not isinstance(target_geometry, gpd.GeoDataFrame):
            raise ValueError("invalid target geometry")

        matches = self.sentinel_polygons[
            self.sentinel_polygons.intersects(target_geometry.to_crs(self.sentinel_polygons.crs).unary_union)]
        matches.rename(columns={"Name": "tile"}, inplace=True)
        tiles = matches[["tile", "geometry"]]

        if calculate_area or eliminate_redundancy:
            centroid = target_geometry.to_crs("EPSG:4326").unary_union.centroid
            lon = centroid.x
            lat = centroid.y
            projection = UTM_proj4_from_latlon(lat, lon)
            tiles_UTM = tiles.to_crs(projection)
            target_UTM = target_geometry.to_crs(projection)
            tiles_UTM["area"] = gpd.overlay(tiles_UTM, target_UTM).geometry.area
            # overlap = gpd.overlay(tiles_UTM, target_UTM)
            # area = overlap.geometry.area

            if eliminate_redundancy:
                # tiles_UTM["area"] = np.array(area)
                tiles_UTM.sort_values(by="area", ascending=False, inplace=True)
                tiles_UTM.reset_index(inplace=True)
                tiles_UTM = tiles_UTM[["tile", "area", "geometry"]]
                remaining_target = target_UTM.unary_union
                remaining_target_area = remaining_target.area
                indices = []

                for i, (tile, area, geometry) in tiles_UTM.iterrows():
                    remaining_target = remaining_target - geometry
                    previous_area = remaining_target_area
                    remaining_target_area = remaining_target.area
                    change_in_area = remaining_target_area - previous_area

                    if change_in_area != 0:
                        indices.append(i)

                    if remaining_target_area == 0:
                        break

                tiles_UTM = tiles_UTM.iloc[indices, :]
                tiles = tiles_UTM.to_crs(tiles.crs)
                tiles.sort_values(by="tile", ascending=True, inplace=True)
                tiles = tiles[["tile", "area", "geometry"]]
            else:
                # tiles["area"] = np.array(area)
                tiles = tiles[["tile", "area", "geometry"]]

        return tiles

    def grid(self, tile: str, cell_size: float = None, buffer=0) -> RasterGrid:
        if cell_size is None:
            cell_size = self.target_resolution

        bbox = self.bbox(tile).buffer(buffer)
        projection = self.UTM_proj4(tile)
        grid = RasterGrid.from_bbox(bbox=bbox, cell_size=cell_size, crs=projection)
        # logger.info(f"tile {cl.place(tile)} at resolution {cl.val(cell_size_degrees)} {grid.shape}")

        return grid

    def land(self, tile: str) -> bool:
        return self.sentinel_polygons[self.sentinel_polygons["Name"].apply(lambda name: name == tile)]["Land"].iloc[0]

    def centroid(self, tile: str) -> shapely.geometry.Point:
        return self.footprint(tile).centroid

    def tile_grids(
            self,
            target_geometry: shapely.geometry.shape or gpd.GeoDataFrame,
            eliminate_redundancy: bool = True) -> gpd.GeoDataFrame:
        tiles = self.tile_footprints(
            target_geometry=target_geometry,
            eliminate_redundancy=eliminate_redundancy,
        )

        tiles["grid"] = tiles["tile"].apply(lambda tile: self.grid(tile))
        tiles = tiles[["tile", "area", "grid", "geometry"]]

        return tiles

sentinel_tile_grid = SentinelTileGrid()
