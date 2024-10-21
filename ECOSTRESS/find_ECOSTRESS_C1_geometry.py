import logging
import sys
from datetime import date, datetime, timedelta
from math import floor
from os.path import exists, join, abspath, dirname
from typing import Union

import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import shapely
import shapely.wkt
from dateutil import parser
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry

import colored_logging as cl
from rasters import RasterGeometry, VectorGeometry
from sentinel_tile_grid import SentinelTileGrid
from ECOSTRESS.find_ECOSTRESS_C1_scene import find_ECOSTRESS_C1_scene


BBOXES_FILENAMES = [
    "/home/gerardo/ecostress/eco_lat_long_info/eco_lat_long_scene_info.csv",
    join(abspath(dirname(__file__)), "eco_lat_long_scene_info.csv")
]

logger = logging.getLogger(__name__)

def evaluate_value(value):
    try:
        return int(value)
    except:
        pass

    try:
        return float(value)
    except:
        pass

    if value in ('True', 'true', 'TRUE'):
        return True

    if value in ('False', 'false', 'FALSE'):
        return False

    try:
        return value.strip().strip('"')
    except:
        pass

    return value


def parse_field(field):
    if isinstance(field, bytes):
        value = field.decode()
    else:
        value = str(field)

    value = evaluate_value(value)

    return value


def UTC_to_solar(time_UTC: datetime, lon: float) -> datetime:
    return time_UTC + timedelta(hours=(np.radians(lon) / np.pi * 12))


def UTM_proj4_from_latlon(lat: float, lon: float) -> str:
    UTM_zone = (floor((lon + 180) / 6) % 60) + 1
    UTM_proj4 = f"+proj=utm +zone={UTM_zone} {'+south ' if lat < 0 else ''}+ellps=WGS84 +datum=WGS84 +units=m +no_defs"

    return UTM_proj4


def load_ECOSTRESS_C1_bboxes(
        start_date_UTC: Union[date, str] = None,
        end_date_UTC: Union[date, str] = None,
        start_date_solar: Union[date, str] = None,
        end_date_solar: Union[date, str] = None,
        hour_min: float = None,
        hour_max: float = None,
        target_geometry_latlon: BaseGeometry = None,
        bboxes_filenames: str = BBOXES_FILENAMES):
    if isinstance(target_geometry_latlon, VectorGeometry):
        target_geometry_latlon = target_geometry_latlon.geometry

    df = None

    for filename in bboxes_filenames:
        if exists(filename):
            logger.info(f"reading ECOSTRESS bounding boxes: {filename}")
            df = pd.read_csv(filename)
            break
        else:
            logger.warning(f"ECOSTRESS bounding boxes file not found: {filename}")

    df.columns = ["bbox_ID", "timestamp", "orbit", "scene", "nw_lat", "nw_lon", "ne_lat", "ne_lon", "sw_lat", "sw_lon",
                  "se_lat", "se_lon"]
    df["time_UTC"] = df.timestamp.apply(lambda timestamp: parser.parse(timestamp))
    geometry = df.apply(lambda row: Polygon(
        [(row.nw_lon, row.nw_lat), (row.ne_lon, row.ne_lat), (row.se_lon, row.se_lat), (row.sw_lon, row.sw_lat)]),
                        axis=1)
    gdf = gpd.GeoDataFrame(df[["time_UTC", "orbit", "scene"]], geometry=geometry, crs="EPSG:4326")
    gdf["time_solar"] = gdf.apply(lambda row: UTC_to_solar(row.time_UTC, row.geometry.centroid.x), axis=1)
    gdf["hour"] = gdf.time_solar.apply(
        lambda time_solar: time_solar.hour + time_solar.minute / 60 + time_solar.second / 3600)
    gdf = gdf[["time_UTC", "time_solar", "hour", "orbit", "scene", "geometry"]]
    gdf = gdf.sort_values(by="time_solar")

    if start_date_UTC is not None:
        if isinstance(start_date_UTC, str):
            start_date_UTC = parser.parse(start_date_UTC).date()

        gdf = gdf[gdf.time_UTC.apply(lambda time_UTC: time_UTC.date() >= start_date_UTC)]

    if end_date_UTC is not None:
        if isinstance(end_date_UTC, str):
            end_date_UTC = parser.parse(end_date_UTC).date()

        gdf = gdf[gdf.time_UTC.apply(lambda time_UTC: time_UTC.date() <= end_date_UTC)]

    if start_date_solar is not None:
        if isinstance(start_date_solar, str):
            start_date_solar = parser.parse(start_date_solar).date()

        gdf = gdf[gdf.time_solar.apply(lambda time_solar: time_solar.date() >= start_date_solar)]

    if end_date_solar is not None:
        if isinstance(end_date_solar, str):
            end_date_solar = parser.parse(end_date_solar).date()

        gdf = gdf[gdf.time_solar.apply(lambda time_solar: time_solar.date() <= end_date_solar)]

    if hour_min is not None:
        gdf = gdf[gdf.hour >= hour_min]

    if hour_max is not None:
        gdf = gdf[gdf.hour <= hour_max]

    if "geometry" not in gdf.columns:
        raise ValueError(f"geometry not in ECOSTRESS scene table columns: {', '.join(gdf.columns)}")

    if target_geometry_latlon is not None:
        gdf = gdf[gdf.geometry.intersects(target_geometry_latlon)]

    return gdf


def find_ECOSTRESS_C1_geometry(
        target_geometry: Union[BaseGeometry, str],
        start_date_UTC: Union[date, str] = None,
        end_date_UTC: Union[date, str] = None,
        start_date_solar: Union[date, str] = None,
        end_date_solar: Union[date, str] = None,
        hour_min: float = None,
        hour_max: float = None,
        max_cloud_percent: float = None,
        results_filename: str = None,
        bboxes_filenames: str = None) -> gpd.GeoDataFrame:
    logger = logging.getLogger(__name__)

    if bboxes_filenames is None:
        bboxes_filenames = BBOXES_FILENAMES

    if isinstance(target_geometry, str):
        if target_geometry.startswith("sentinel:"):
            tile = target_geometry.split(":")[-1]
            sentinel_tile_grid = SentinelTileGrid()
            target_geometry = sentinel_tile_grid.footprint(tile).geometry
        elif exists(target_geometry):
            target_geometry = gpd.read_file(target_geometry).to_crs("EPSG:4326").unary_union
        else:
            target_geometry = shapely.wkt.loads(target_geometry)
    elif isinstance(target_geometry, RasterGeometry):
        target_geometry = target_geometry.corner_polygon_latlon.geometry

    logger.info("searching ECOSTRESS C1 bounding boxes")

    gdf = load_ECOSTRESS_C1_bboxes(
        start_date_UTC=start_date_UTC,
        end_date_UTC=end_date_UTC,
        start_date_solar=start_date_solar,
        end_date_solar=end_date_solar,
        hour_min=hour_min,
        hour_max=hour_max,
        target_geometry_latlon=target_geometry,
        bboxes_filenames=bboxes_filenames
    )

    logger.info(f"found {cl.val(len(gdf))} matching bounding boxes")

    matching_rows = []
    matching_geometry = []

    for i, row in gdf.iterrows():
        time_UTC = row.time_UTC
        date_UTC = row.time_UTC.date()
        time_solar = row.time_solar
        hour = row.hour
        orbit = row.orbit
        scene = row.scene
        bbox_polygon = row.geometry

        logger.info(
            "time: " + cl.time(f"{time_UTC:%Y-%m-%d %H:%M:%S} UTC") +
            " hour: " + cl.time(f"{hour:0.2f}") +
            " orbit: " + cl.val(orbit) +
            " scene: " + cl.val(scene)
        )

        logger.info("searching files")
        try:
            filenames_dict = find_ECOSTRESS_C1_scene(
                orbit=orbit,
                scene=scene,
                date_UTC=date_UTC
            )
        except Exception as e:
            logger.warning(e)
            continue

        if "L1B_GEO" not in filenames_dict:
            logger.warning(f"cannot find L1B GEO file for orbit {orbit} scene {scene}")
            continue

        L1B_GEO_filename = filenames_dict["L1B_GEO"]

        if "L1B_RAD" not in filenames_dict:
            logger.warning(f"cannot find L1B RAD file for orbit {orbit} scene {scene}")
            continue

        L1B_RAD_filename = filenames_dict["L1B_RAD"]

        if "L2_LSTE" not in filenames_dict:
            logger.warning(f"cannot find L2 LSTE file for orbit {orbit} scene {scene}")
            continue

        L2_LSTE_filename = filenames_dict["L2_LSTE"]

        if "L2_CLOUD" not in filenames_dict:
            logger.warning(f"cannot find L2 CLOUD file for orbit {orbit} scene {scene}")
            continue

        L2_CLOUD_filename = filenames_dict["L2_CLOUD"]

        with h5py.File(L1B_GEO_filename, "r") as file:
            day = parse_field(file['StandardMetadata/DayNightFlag'][()]) == 'Day'

            if not day:
                logger.warning(f"skipping non-day orbit: {orbit} scene: {scene} hour: {hour:0.2f}")
                continue

            lon = file["Geolocation/longitude"]
            lat = file["Geolocation/latitude"]

            corner_polygon = Polygon([
                (lon[0, 0], lat[0, 0]),
                (lon[0, lon.shape[1] - 1], lat[0, lat.shape[1] - 1]),
                (lon[lon.shape[0] - 1, lon.shape[1] - 1], lat[lat.shape[0] - 1, lat.shape[1] - 1]),
                (lon[lon.shape[0] - 1, 0], lat[lat.shape[0] - 1, 0])
            ])

        logger.info(f"bbox: {cl.val(shapely.wkt.dumps(bbox_polygon))}")
        logger.info(f"corner polygon: {cl.val(shapely.wkt.dumps(corner_polygon))}")

        if not corner_polygon.intersects(target_geometry):
            logger.warning(f"ECOSTRESS orbit {orbit} scene {scene} does not intersect target geometry")
            continue

        with h5py.File(L2_CLOUD_filename, "r") as file:
            cloud_percent = int(parse_field(file['L2 CLOUD Metadata/QAPercentCloudCover'][()])[1:-1])

        if max_cloud_percent is not None and cloud_percent > max_cloud_percent:
            logger.warning(
                f"ECOSTRESS orbit {orbit} scene {scene} cloud coverage {cloud_percent}% exceeds maximum {max_cloud_percent}%")
            continue

        logger.info(f"cloud: {cl.val(f'{cloud_percent}%')}")

        matching_rows.append([time_UTC, time_solar, hour, orbit, scene, cloud_percent, L1B_GEO_filename, L1B_RAD_filename, L2_LSTE_filename, L2_CLOUD_filename])
        matching_geometry.append(corner_polygon)

    df = pd.DataFrame(matching_rows, columns=["time_UTC", "time_solar", "hour", "orbit", "scene", "cloud_percent", "L1B_GEO", "L1B_RAD", "L2_LSTE", "L2_CLOUD"])
    gdf = gpd.GeoDataFrame(df, geometry=matching_geometry, crs="EPSG:4326")

    target_gdf_latlon = gpd.GeoDataFrame({}, geometry=[target_geometry], crs="EPSG:4326")
    centroid = target_gdf_latlon.unary_union.centroid
    lon = centroid.x
    lat = centroid.y
    projection = UTM_proj4_from_latlon(lat, lon)
    scenes_UTM = gdf.to_crs(projection)
    target_UTM = target_gdf_latlon.to_crs(projection)
    overlap = gpd.overlay(scenes_UTM, target_UTM)
    area = overlap.geometry.area
    gdf["area"] = area
    gdf = gdf[["time_UTC", "time_solar", "hour", "orbit", "scene", "cloud_percent", "L1B_GEO", "L1B_RAD", "L2_LSTE", "L2_CLOUD", "area", "geometry"]]
    logger.info(f"found {cl.val(len(gdf))} scenes matching corner polygons")

    for i, row in gdf.iterrows():
        time_UTC = row.time_UTC
        hour = row.hour
        orbit = row.orbit
        scene = row.scene
        cloud_percent = row.cloud_percent
        area = row.area

        logger.info(
            "* time: " + cl.time(f"{time_UTC:%Y-%m-%d} UTC") +
            " hour: " + cl.time(f"{hour:0.2f}") +
            " orbit: " + cl.val(orbit) +
            " scene: " + cl.val(scene) +
            " cloud: " + cl.val(f"{cloud_percent}%") +
            " area: " + cl.val(f"{area:0.2f} m^2")
        )

    if results_filename is not None:
        gdf.to_file(results_filename, driver="GeoJSON")

    return gdf


def main(argv=sys.argv):
    geometry = argv[1]

    if "--start-date-UTC" in argv:
        start_date_UTC = parser.parse(argv[argv.index("--start-date-UTC") + 1]).date()
    elif "--start" in argv:
        start_date_UTC = parser.parse(argv[argv.index("--start") + 1]).date()
    else:
        start_date_UTC = None

    if "--end-date-UTC" in argv:
        end_date_UTC = parser.parse(argv[argv.index("--end-date-UTC") + 1]).date()
    elif "--end" in argv:
        end_date_UTC = parser.parse(argv[argv.index("--end") + 1]).date()
    else:
        end_date_UTC = None

    if "--max-cloud" in argv:
        max_cloud_percent = float(argv[argv.index("--max-cloud") + 1])
    else:
        max_cloud_percent = None

    find_ECOSTRESS_C1_geometry(
        target_geometry=geometry,
        start_date_UTC=start_date_UTC,
        end_date_UTC=end_date_UTC,
        max_cloud_percent=max_cloud_percent
    )


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
