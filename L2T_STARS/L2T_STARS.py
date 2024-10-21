import logging
import shutil
import socket
import subprocess
import sys
import urllib
from datetime import datetime, timedelta, date
from glob import glob
from os import makedirs, system, remove
from os.path import join, abspath, dirname, expanduser, exists, basename
from shutil import which
from traceback import format_exception
from typing import Union
from uuid import uuid4

import numpy as np
import pandas as pd
from dateutil import parser
from dateutil.rrule import rrule, DAILY
from scipy import stats

import ECOSTRESS
import colored_logging as cl
from ecostress_cmr import CMRServerUnreachable
import rasters as rt
from ECOSTRESS.L2_LSTE import L2TLSTE
from ECOSTRESS.L2_STARS import L2TSTARS
from ECOSTRESS.exit_codes import ANCILLARY_LATENCY, LAND_FILTER, SUCCESS_EXIT_CODE, ECOSTRESSExitCodeException, \
    UnableToParseRunConfig, \
    MissingRunConfigValue, RUNCONFIG_FILENAME_NOT_SUPPLIED, AncillaryServerUnreachable, AncillaryLatency, LandFilter, \
    InputFilesInaccessible, DownloadFailed, ANCILLARY_SERVER_UNREACHABLE, BlankOutput, DOWNLOAD_FAILED, \
    UNCLASSIFIED_FAILURE_EXIT_CODE
from ECOSTRESS.runconfig import ECOSTRESSRunConfig
from ECOSTRESS_colors import NDVI_COLORMAP, ALBEDO_COLORMAP
from HLS import HLSLandsatMissing, HLSSentinelMissing, HLS
from HLS import HLSTileNotAvailable, HLSSentinelNotAvailable, HLSLandsatNotAvailable, HLSDownloadFailed, HLSNotAvailable
from HLS.HLS2 import HLSBandNotAcquired, HLS2CMR, CMR_SEARCH_URL
from LPDAAC.LPDAACDataPool import LPDAACServerUnreachable
from MODLAND import find_MODLAND_tiles
from VIIRS import VIIRSDownloaderAlbedo, VIIRSDownloaderNDVI
from VIIRS.VNP09GA import VNP09GA, VIIRSUnavailableError
from VIIRS.VNP43IA4 import VNP43IA4
from VIIRS.VNP43MA3 import VNP43MA3
from VNP43NRT import VNP43NRT
from daterange import get_date
from rasters import Raster, RasterGeometry
from timer import Timer

with open(join(abspath(dirname(__file__)), "version.txt")) as f:
    version = f.read()

__version__ = version

L2T_STARS_TEMPLATE = join(abspath(dirname(__file__)), "L2T_STARS.xml")
DEFAULT_WORKING_DIRECTORY = "."
DEFAULT_BUILD = "0700"
PRIMARY_VARIABLE = "NDVI"
DEFAULT_OUTPUT_DIRECTORY = "L2T_STARS_output"
DEFAULT_STARS_SOURCES_DIRECTORY = "L2T_STARS_SOURCES"
DEFAULT_STARS_INDICES_DIRECTORY = "L2T_STARS_INDICES"
DEFAULT_STARS_MODEL_DIRECTORY = "L2T_STARS_MODEL"
DEFAULT_STARS_PRODUCTS_DIRECTORY = "STARS_products"
DEFAULT_HLS_DOWNLOAD_DIRECTORY = "HLS2_download"
DEFAULT_LANDSAT_DOWNLOAD_DIRECTORY = "HLS2_download"
DEFAULT_HLS_PRODUCTS_DIRECTORY = "HLS2_products"
DEFAULT_VIIRS_DOWNLOAD_DIRECTORY = "VIIRS_download"
DEFAULT_VIIRS_PRODUCTS_DIRECTORY = "VIIRS_products"
DEFAUL_VIIRS_MOSAIC_DIRECTORY = "VIIRS_mosaic"
DEFAULT_GEOS5FP_DOWNLOAD_DIRECTORY = "GEOS5FP_download"
DEFAULT_GEOS5FP_PRODUCTS_DIRECTORY = "GEOS5FP_products"
DEFAULT_VNP09GA_PRODUCTS_DIRECTORY = "VNP09GA_products"
DEFAULT_VNP43NRT_PRODUCTS_DIRECTORY = "VNP43NRT_products"
VIIRS_GIVEUP_DAYS = 4
DEFAULT_SPINUP_DAYS = 7
DEFAULT_TARGET_RESOLUTION = 70
DEFAULT_NDVI_RESOLUTION = 490
DEFAULT_ALBEDO_RESOLUTION = 980
DEFAULT_USE_SPATIAL = False
DEFAULT_USE_VNP43NRT = True
DEFAULT_CALIBRATE_FINE = False

L2T_STARS_SHORT_NAME = "ECO_L2T_STARS"
L2T_STARS_LONG_NAME = "ECOSTRESS Tiled Ancillary NDVI and Albedo L2 Global 70 m"

logger = logging.getLogger(__name__)

class Prior:
    def __init__(
            self,
            using_prior: bool = False,
            prior_date_UTC: date = None,
            L2T_STARS_prior_filename: str = None,
            prior_NDVI_filename: str = None,
            prior_NDVI_UQ_filename: str = None,
            prior_NDVI_bias_filename: str = None,
            prior_NDVI_bias_UQ_filename: str = None,
            prior_albedo_filename: str = None,
            prior_albedo_UQ_filename: str = None,
            prior_albedo_bias_filename: str = None,
            prior_albedo_bias_UQ_filename: str = None):
        self.using_prior = using_prior
        self.prior_date_UTC = prior_date_UTC
        self.L2T_STARS_prior_filename = L2T_STARS_prior_filename
        self.prior_NDVI_filename = prior_NDVI_filename
        self.prior_NDVI_UQ_filename = prior_NDVI_UQ_filename
        self.prior_NDVI_bias_filename = prior_NDVI_bias_filename
        self.prior_NDVI_bias_UQ_filename = prior_NDVI_bias_UQ_filename
        self.prior_albedo_filename = prior_albedo_filename
        self.prior_albedo_UQ_filename = prior_albedo_UQ_filename
        self.prior_albedo_bias_filename = prior_albedo_bias_filename
        self.prior_albedo_bias_UQ_filename = prior_albedo_bias_UQ_filename


def generate_L2T_STARS_runconfig(
        L2T_LSTE_filename: str,
        prior_L2T_STARS_filename: str = "",
        orbit: int = None,
        scene: int = None,
        tile: str = None,
        time_UTC: Union[datetime, str] = None,
        working_directory: str = None,
        sources_directory: str = None,
        indices_directory: str = None,
        model_directory: str = None,
        executable_filename: str = None,
        output_directory: str = None,
        runconfig_filename: str = None,
        log_filename: str = None,
        build: str = None,
        processing_node: str = None,
        production_datetime: datetime = None,
        job_ID: str = None,
        instance_ID: str = None,
        product_counter: int = None,
        template_filename: str = None) -> str:
    logger.info(
        f"started generating L2T_STARS run-config for orbit {cl.val(orbit)} scene {cl.val(scene)}")
    timer = Timer()

    L2T_LSTE_granule = L2TLSTE(L2T_LSTE_filename)

    if orbit is None:
        orbit = L2T_LSTE_granule.orbit

    if scene is None:
        scene = L2T_LSTE_granule.scene

    if tile is None:
        tile = L2T_LSTE_granule.tile

    if time_UTC is None:
        time_UTC = L2T_LSTE_granule.time_UTC

    if build is None:
        build = DEFAULT_BUILD

    pattern = join(working_directory,
                   f"ECOv002_L2T_STARS_{orbit:05d}_{scene:03d}_{tile}_*_{build}_*.xml")
    logger.info(f"scanning for previous run-configs: {cl.val(pattern)}")
    previous_runconfigs = glob(pattern)
    previous_runconfig_count = len(previous_runconfigs)

    if previous_runconfig_count > 0:
        logger.info(
            f"found {cl.val(previous_runconfig_count)} previous run-configs")
        previous_runconfig = sorted(previous_runconfigs)[-1]
        logger.info(f"previous run-config: {cl.file(previous_runconfig)}")

        return previous_runconfig

    if template_filename is None:
        template_filename = L2T_STARS_TEMPLATE

    template_filename = abspath(expanduser(template_filename))

    if production_datetime is None:
        production_datetime = datetime.utcnow()

    if product_counter is None:
        product_counter = 1

    timestamp = f"{time_UTC:%Y%m%d}"
    granule_ID = f"ECOv002_L2T_STARS_{orbit:05d}_{scene:03d}_{tile}_{timestamp}_{build}_{product_counter:02d}"

    if runconfig_filename is None:
        runconfig_filename = join(
            working_directory, "runconfig", f"{granule_ID}.xml")

    runconfig_filename = abspath(expanduser(runconfig_filename))

    if exists(runconfig_filename):
        logger.info(f"run-config already exists {cl.file(runconfig_filename)}")
        return runconfig_filename

    if working_directory is None:
        working_directory = granule_ID

    working_directory = abspath(expanduser(working_directory))

    if sources_directory is None:
        sources_directory = join(
            working_directory, DEFAULT_STARS_SOURCES_DIRECTORY)

    if indices_directory is None:
        indices_directory = join(
            working_directory, DEFAULT_STARS_INDICES_DIRECTORY)

    if model_directory is None:
        model_directory = join(
            working_directory, DEFAULT_STARS_MODEL_DIRECTORY)

    if executable_filename is None:
        executable_filename = which("L2T_STARS")

    if executable_filename is None:
        executable_filename = "L2T_STARS"

    if output_directory is None:
        output_directory = join(working_directory, DEFAULT_OUTPUT_DIRECTORY)

    output_directory = abspath(expanduser(output_directory))

    if log_filename is None:
        log_filename = join(working_directory, "log", f"{granule_ID}.log")

    log_filename = abspath(expanduser(log_filename))

    if processing_node is None:
        processing_node = socket.gethostname()

    if job_ID is None:
        job_ID = timestamp

    if instance_ID is None:
        instance_ID = str(uuid4())

    L2T_LSTE_filename = abspath(expanduser(L2T_LSTE_filename))

    logger.info(f"loading L2T_STARS template: {cl.file(template_filename)}")

    with open(template_filename, "r") as file:
        template = file.read()

    logger.info(f"orbit: {cl.val(orbit)}")
    template = template.replace("orbit_number", f"{orbit:05d}")
    logger.info(f"scene: {cl.val(scene)}")
    template = template.replace("scene_ID", f"{scene:03d}")
    logger.info(f"tile: {cl.val(tile)}")
    template = template.replace("tile_ID", f"{tile}")
    logger.info(f"L2T_LSTE file: {cl.file(L2T_LSTE_filename)}")
    template = template.replace("L2T_LSTE_filename", L2T_LSTE_filename)
    logger.info(f"prior L2T_STARS file: {cl.file(prior_L2T_STARS_filename)}")
    template = template.replace("prior_L2T_STARS_filename", prior_L2T_STARS_filename)
    logger.info(f"working directory: {cl.dir(working_directory)}")
    template = template.replace("working_directory", working_directory)
    logger.info(f"sources directory: {cl.dir(sources_directory)}")
    template = template.replace("sources_directory", sources_directory)
    logger.info(f"indices directory: {cl.dir(indices_directory)}")
    template = template.replace("indices_directory", indices_directory)
    logger.info(f"model directory: {cl.dir(model_directory)}")
    template = template.replace("model_directory", model_directory)
    logger.info(f"executable: {cl.file(executable_filename)}")
    template = template.replace("executable_filename", executable_filename)
    logger.info(f"output directory: {cl.dir(output_directory)}")
    template = template.replace("output_directory", output_directory)
    logger.info(f"run-config: {cl.file(runconfig_filename)}")
    template = template.replace("runconfig_filename", runconfig_filename)
    logger.info(f"log: {cl.file(log_filename)}")
    template = template.replace("log_filename", log_filename)
    logger.info(f"build: {cl.val(build)}")
    template = template.replace("build_ID", build)
    logger.info(f"processing node: {cl.val(processing_node)}")
    template = template.replace("processing_node", processing_node)
    logger.info(f"production date/time: {cl.time(production_datetime)}")
    template = template.replace("production_datetime", timestamp)
    logger.info(f"job ID: {cl.val(job_ID)}")
    template = template.replace("job_ID", job_ID)
    logger.info(f"instance ID: {cl.val(instance_ID)}")
    template = template.replace("instance_ID", instance_ID)
    logger.info(f"product counter: {cl.val(product_counter)}")
    template = template.replace("product_counter", f"{product_counter:02d}")

    makedirs(dirname(abspath(runconfig_filename)), exist_ok=True)
    logger.info(f"writing run-config file: {cl.file(runconfig_filename)}")

    with open(runconfig_filename, "w") as file:
        file.write(template)

    logger.info(
        f"finished generating L2T_STARS run-config for orbit {cl.val(orbit)} scene {cl.val(scene)} ({timer})")

    return runconfig_filename


class L2TSTARSConfig(ECOSTRESSRunConfig):
    def __init__(self, filename: str):
        logger.info(f"loading L2T_STARS run-config: {cl.file(filename)}")
        runconfig = self.read_runconfig(filename)

        # print(JSON_highlight(runconfig))

        try:
            if "StaticAncillaryFileGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup in L2T_STARS run-config: {filename}")

            if "L2T_STARS_WORKING" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2T_STARS_WORKING in L2T_STARS run-config: {filename}")

            working_directory = abspath(
                runconfig["StaticAncillaryFileGroup"]["L2T_STARS_WORKING"])
            logger.info(f"working directory: {cl.dir(working_directory)}")

            if "L2T_STARS_SOURCES" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2T_STARS_WORKING in L2T_STARS run-config: {filename}")

            sources_directory = abspath(
                runconfig["StaticAncillaryFileGroup"]["L2T_STARS_SOURCES"])
            logger.info(f"sources directory: {cl.dir(sources_directory)}")

            if "L2T_STARS_INDICES" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2T_STARS_INDICES in L2T_STARS run-config: {filename}")

            indices_directory = abspath(
                runconfig["StaticAncillaryFileGroup"]["L2T_STARS_INDICES"])
            logger.info(f"indices directory: {cl.dir(indices_directory)}")

            if "L2T_STARS_MODEL" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2T_STARS_MODEL in L2T_STARS run-config: {filename}")

            model_directory = abspath(
                runconfig["StaticAncillaryFileGroup"]["L2T_STARS_MODEL"])
            logger.info(f"model directory: {cl.dir(model_directory)}")

            if "ProductPathGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup in L2T_STARS run-config: {filename}")

            if "ProductPath" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductPath in L2T_STARS run-config: {filename}")

            output_directory = abspath(
                runconfig["ProductPathGroup"]["ProductPath"])
            logger.info(f"output directory: {cl.dir(output_directory)}")

            if "InputFileGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup in L2G_L2T_LSTE run-config: {filename}")

            if "L2T_LSTE" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L2T_LSTE in L2T_STARS run-config: {filename}")

            L2T_LSTE_filename = abspath(
                runconfig["InputFileGroup"]["L2T_LSTE"])
            logger.info(f"L2T_LSTE file: {cl.file(L2T_LSTE_filename)}")

            if "L2T_STARS_PRIOR" in runconfig["InputFileGroup"]:
                L2T_STARS_prior_filename = runconfig["InputFileGroup"]["L2T_STARS_PRIOR"]

                if L2T_STARS_prior_filename != "" and exists(L2T_STARS_prior_filename):
                    L2T_STARS_prior_filename = abspath(L2T_STARS_prior_filename)

                logger.info(f"L2T_STARS prior file: {cl.file(L2T_STARS_prior_filename)}")
            else:
                L2T_STARS_prior_filename = None

            orbit = int(runconfig["Geometry"]["OrbitNumber"])
            logger.info(f"orbit: {cl.val(orbit)}")

            if "SceneId" not in runconfig["Geometry"]:
                raise MissingRunConfigValue(
                    f"missing Geometry/SceneId in L2T_STARS run-config: {filename}")

            scene = int(runconfig["Geometry"]["SceneId"])
            logger.info(f"scene: {cl.val(scene)}")

            if "TileId" not in runconfig["Geometry"]:
                raise MissingRunConfigValue(
                    f"missing Geometry/TileId in L2T_STARS run-config: {filename}")

            tile = str(runconfig["Geometry"]["TileId"])
            logger.info(f"tile: {cl.val(tile)}")

            if "ProductionDateTime" not in runconfig["JobIdentification"]:
                raise MissingRunConfigValue(
                    f"missing JobIdentification/ProductionDateTime in L2T_STARS run-config {filename}")

            production_datetime = parser.parse(
                runconfig["JobIdentification"]["ProductionDateTime"])
            logger.info(f"production time: {cl.time(production_datetime)}")

            if "BuildID" not in runconfig["PrimaryExecutable"]:
                raise MissingRunConfigValue(
                    f"missing PrimaryExecutable/BuildID in L2T_STARS run-config {filename}")

            build = str(runconfig["PrimaryExecutable"]["BuildID"])

            if "ProductCounter" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductCounter in L2T_STARS run-config {filename}")

            product_counter = int(
                runconfig["ProductPathGroup"]["ProductCounter"])

            L2T_LSTE_granule = L2TLSTE(L2T_LSTE_filename)
            time_UTC = L2T_LSTE_granule.time_UTC

            granule_ID = f"ECOv002_L2T_STARS_{tile}_{time_UTC:%Y%m%d}_{build}_{product_counter:02d}"

            L2T_STARS_granule_directory = join(output_directory, granule_ID)
            L2T_STARS_zip_filename = f"{L2T_STARS_granule_directory}.zip"
            L2T_STARS_browse_filename = f"{L2T_STARS_granule_directory}.png"

            self.working_directory = working_directory
            self.sources_directory = sources_directory
            self.indices_directory = indices_directory
            self.model_directory = model_directory
            self.output_directory = output_directory
            self.L2T_LSTE_filename = L2T_LSTE_filename
            self.L2T_STARS_prior_filename = L2T_STARS_prior_filename
            self.orbit = orbit
            self.scene = scene
            self.tile = tile
            self.production_datetime = production_datetime
            self.build = build
            self.product_counter = product_counter
            self.granule_ID = granule_ID
            self.L2T_STARS_granule_directory = L2T_STARS_granule_directory
            self.L2T_STARS_zip_filename = L2T_STARS_zip_filename
            self.L2T_STARS_browse_filename = L2T_STARS_browse_filename

        except MissingRunConfigValue as e:
            raise e
        except ECOSTRESSExitCodeException as e:
            raise e
        except Exception as e:
            logger.exception(e)
            raise UnableToParseRunConfig(
                f"unable to parse run-config file: {filename}")

def generate_filename(directory: str, variable: str, date_UTC: Union[date, str], tile: str, cell_size: int) -> str:
    if isinstance(date_UTC, str):
        date_UTC = parser.parse(date_UTC).date()

    variable = str(variable)
    timestamp = date_UTC.strftime("%Y-%m-%d")
    tile = str(tile)
    cell_size = int(cell_size)
    filename = join(directory, f"STARS_{variable}_{timestamp}_{tile}_{cell_size}m.tif")
    makedirs(dirname(filename), exist_ok=True)

    return filename

def calibrate_fine_to_coarse(fine_image: Raster, coarse_image: Raster) -> Raster:
    aggregated_image = fine_image.to_geometry(coarse_image.geometry, resampling="average")
    x = np.array(aggregated_image).flatten()
    y = np.array(coarse_image).flatten()
    mask = ~np.isnan(x) & ~np.isnan(y)

    if np.count_nonzero(mask) < 30:
        return fine_image

    x = x[mask]
    y = y[mask]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    calibrated_image = fine_image * slope + intercept

    return calibrated_image

def generate_NDVI_coarse_image(date_UTC: Union[date, str], VIIRS_connection: VIIRSDownloaderNDVI, geometry: RasterGeometry = None) -> Raster:
    coarse_image = VIIRS_connection.NDVI(date_UTC=date_UTC, geometry=geometry)
    coarse_image = rt.where(coarse_image == 0, np.nan, coarse_image)

    return coarse_image

def generate_NDVI_fine_image(date_UTC: Union[date, str], tile: str, HLS_connection: HLS) -> Raster:
    fine_image = HLS_connection.NDVI(tile=tile, date_UTC=date_UTC)
    fine_image = rt.where(fine_image == 0, np.nan, fine_image)

    return fine_image

def generate_albedo_coarse_image(date_UTC: Union[date, str], VIIRS_connection: VIIRSDownloaderAlbedo, geometry: RasterGeometry = None) -> Raster:
    coarse_image = VIIRS_connection.albedo(date_UTC=date_UTC, geometry=geometry)
    coarse_image = rt.where(coarse_image == 0, np.nan, coarse_image)

    return coarse_image

def generate_albedo_fine_image(date_UTC: Union[date, str], tile: str, HLS_connection: HLS) -> Raster:
    fine_image = HLS_connection.albedo(tile=tile, date_UTC=date_UTC)
    fine_image = rt.where(fine_image == 0, np.nan, fine_image)

    return fine_image

def generate_input_staging_directory(input_staging_directory: str, tile: str, prefix: str) -> str:
    directory = join(input_staging_directory, f"{prefix}_{tile}")
    makedirs(directory, exist_ok=True)

    return directory

def generate_NDVI_coarse_directory(input_staging_directory: str, tile: str) -> str:
    return generate_input_staging_directory(input_staging_directory, tile, "NDVI_coarse")

def generate_NDVI_fine_directory(input_staging_directory: str, tile: str) -> str:
    return generate_input_staging_directory(input_staging_directory, tile, "NDVI_fine")

def generate_albedo_coarse_directory(input_staging_directory: str, tile: str) -> str:
    return generate_input_staging_directory(input_staging_directory, tile, "albedo_coarse")

def generate_albedo_fine_directory(input_staging_directory: str, tile: str) -> str:
    return generate_input_staging_directory(input_staging_directory, tile, "albedo_fine")

def generate_output_directory(working_directory: str, date_UTC: Union[date, str], tile: str) -> str:
    if isinstance(date_UTC, str):
        date_UTC = parser.parse(date_UTC).date()

    directory = join(working_directory, f"julia_output_{date_UTC:%y.%m.%d}_{tile}")
    makedirs(directory, exist_ok=True)

    return directory

def generate_model_state_tile_date_directory(model_directory: str, tile: str, date_UTC: Union[date, str]) -> str:
    if isinstance(date_UTC, str):
        date_UTC = parser.parse(date_UTC).date()

    directory = join(model_directory, tile, f"{date_UTC:%Y-%m-%d}")
    makedirs(directory, exist_ok=True)

    return directory

# def process_julia_data_fusion_no_prior(
#         tile: str,
#         coarse_cell_size: int,
#         fine_cell_size: int,
#         VIIRS_start_date: date,
#         VIIRS_end_date: date,
#         HLS_start_date: date,
#         HLS_end_date: date,
#         coarse_directory: str,
#         fine_directory: str,
#         posterior_filename: str,
#         posterior_UQ_filename: str,
#         posterior_bias_filename: str,
#         posterior_bias_UQ_filename: str):
#     julia_script_filename = join(abspath(dirname(__file__)), "process_ECOSTRESS_data_fusion.jl")
#     STARS_source_directory = join(abspath(dirname(__file__)), "..", "STARS_jl")
#     command = f'cd "{STARS_source_directory}" && julia --project=@. "{julia_script_filename}" "{tile}" "{coarse_cell_size}" "{fine_cell_size}" "{VIIRS_start_date}" "{VIIRS_end_date}" "{HLS_start_date}" "{HLS_end_date}" "{coarse_directory}" "{fine_directory}" "{posterior_filename}" "{posterior_UQ_filename}" "{posterior_bias_filename}" "{posterior_bias_UQ_filename}"'
#     logger.info(command)
#     system(command)

def process_julia_data_fusion(
        tile: str,
        coarse_cell_size: int,
        fine_cell_size: int,
        VIIRS_start_date: date,
        VIIRS_end_date: date,
        HLS_start_date: date,
        HLS_end_date: date,
        coarse_directory: str,
        fine_directory: str,
        posterior_filename: str,
        posterior_UQ_filename: str,
        posterior_bias_filename: str,
        posterior_bias_UQ_filename: str,
        prior_filename: str = None,
        prior_UQ_filename: str = None,
        prior_bias_filename: str = None,
        prior_bias_UQ_filename: str = None):
    julia_script_filename = join(abspath(dirname(__file__)), "process_ECOSTRESS_data_fusion.jl")
    STARS_source_directory = join(abspath(dirname(__file__)), "..", "STARS_jl")
    # command = f'cd "{STARS_source_directory}" && julia --project=@ECOSTRESS "{julia_script_filename}" "{tile}" "{coarse_cell_size}" "{fine_cell_size}" "{VIIRS_start_date}" "{VIIRS_end_date}" "{HLS_start_date}" "{HLS_end_date}" "{coarse_directory}" "{fine_directory}" "{posterior_filename}" "{posterior_UQ_filename}" "{posterior_bias_filename}" "{posterior_bias_UQ_filename}"'
    command = f'export JULIA_NUM_THREADS=auto; julia --threads auto "{julia_script_filename}" "{tile}" "{coarse_cell_size}" "{fine_cell_size}" "{VIIRS_start_date}" "{VIIRS_end_date}" "{HLS_start_date}" "{HLS_end_date}" "{coarse_directory}" "{fine_directory}" "{posterior_filename}" "{posterior_UQ_filename}" "{posterior_bias_filename}" "{posterior_bias_UQ_filename}"'

    if all([filename is not None and exists(filename) for filename in [prior_filename, prior_UQ_filename, prior_bias_filename, prior_bias_UQ_filename]]):
        logger.info("passing prior into Julia data fusion system")
        command += f' "{prior_filename}" "{prior_UQ_filename}" "{prior_bias_filename}" "{prior_bias_UQ_filename}"'

    logger.info(command)
    # system(command)
    subprocess.run(command, shell=True)

def retrieve_STARS_sources(
        tile: str,
        geometry: RasterGeometry,
        HLS_start_date: date,
        HLS_end_date: date,
        VIIRS_start_date: date,
        VIIRS_end_date: date,
        HLS_connection: HLS2CMR,
        VIIRS_connection: VNP43NRT):
    logger.info(
        f"retrieving HLS sources for tile {cl.place(tile)} from {cl.time(HLS_start_date)} to {cl.time(HLS_end_date)}")
    for processing_date in [get_date(dt) for dt in rrule(DAILY, dtstart=HLS_start_date, until=HLS_end_date)]:
        try:
            logger.info(
                f"retrieving HLS Sentinel at tile {cl.place(tile)} on date {cl.time(processing_date)}")
            HLS_connection.sentinel(tile=tile, date_UTC=processing_date)
            logger.info(
                f"retrieving HLS Landsat at tile {cl.place(tile)} on date {cl.time(processing_date)}")
            HLS_connection.landsat(tile=tile, date_UTC=processing_date)
        except HLSDownloadFailed as e:
            logger.exception(e)
            raise DownloadFailed(e)
        except HLSTileNotAvailable as e:
            logger.warning(e)
        except HLSSentinelNotAvailable as e:
            logger.warning(e)
        except HLSLandsatNotAvailable as e:
            logger.warning(e)
        except Exception as e:
            logger.warning("Exception raised while retrieving HLS tiles")
            logger.exception(e)
            continue

    # VIIRS_start_date = start_date - timedelta(days=16)
    logger.info(
        f"retrieving VIIRS sources for tile {cl.place(tile)} from {cl.time(VIIRS_start_date)} to {cl.time(VIIRS_end_date)}")

    VIIRS_connection.prefetch_VNP09GA(
        start_date=VIIRS_start_date,
        end_date=VIIRS_end_date,
        geometry=geometry,
    )

def generate_STARS_inputs(
        tile: str,
        date_UTC: date,
        HLS_start_date: date,
        HLS_end_date: date,
        VIIRS_start_date: date,
        VIIRS_end_date: date,
        NDVI_resolution: int,
        albedo_resolution: int,
        target_resolution: int,
        NDVI_coarse_geometry: RasterGeometry,
        albedo_coarse_geometry: RasterGeometry,
        working_directory: str,
        NDVI_coarse_directory: str,
        NDVI_fine_directory: str,
        albedo_coarse_directory: str,
        albedo_fine_directory: str,
        HLS_connection: HLS2CMR,
        NDVI_VIIRS_connection: VIIRSDownloaderNDVI,
        albedo_VIIRS_connection: VIIRSDownloaderAlbedo,
        calibrate_fine: True):
    missing_coarse_dates = set()

    for processing_date in [get_date(dt) for dt in rrule(DAILY, dtstart=VIIRS_start_date, until=VIIRS_end_date)]:
        logger.info(f"preparing coarse image for STARS NDVI at {cl.place(tile)} on {cl.time(processing_date)}")

        try:
            NDVI_coarse_image = generate_NDVI_coarse_image(
                date_UTC=processing_date,
                VIIRS_connection=NDVI_VIIRS_connection,
                geometry=NDVI_coarse_geometry
            )

            NDVI_coarse_filename = generate_filename(
                directory=NDVI_coarse_directory,
                variable="NDVI",
                date_UTC=processing_date,
                tile=tile,
                cell_size=NDVI_resolution
            )

            logger.info(f"saving coarse image for STARS NDVI at {cl.place(tile)} on {cl.time(processing_date)}: {NDVI_coarse_filename}")
            NDVI_coarse_image.to_geotiff(NDVI_coarse_filename)

            if processing_date >= HLS_start_date:
                logger.info(f"preparing fine image for STARS NDVI at {cl.place(tile)} on {cl.time(processing_date)}")

                try:
                    NDVI_fine_image = generate_NDVI_fine_image(
                        date_UTC=processing_date,
                        tile=tile,
                        HLS_connection=HLS_connection
                    )

                    if calibrate_fine:
                        logger.info(f"calibrating fine image for STARS NDVI at {cl.place(tile)} on {cl.time(processing_date)}")
                        NDVI_fine_image = calibrate_fine_to_coarse(NDVI_fine_image, NDVI_coarse_image)

                    NDVI_fine_filename = generate_filename(
                        directory=NDVI_fine_directory,
                        variable="NDVI",
                        date_UTC=processing_date,
                        tile=tile,
                        cell_size=target_resolution
                    )

                    logger.info(f"saving fine image for STARS NDVI at {cl.place(tile)} on {cl.time(processing_date)}: {NDVI_fine_filename}")
                    NDVI_fine_image.to_geotiff(NDVI_fine_filename)
                except:
                    logger.info(f"HLS is not available on {processing_date}")
        except Exception as e:
            logger.exception(e)
            logger.warning(f"unable to produce coarse NDVI for date {processing_date}")
            missing_coarse_dates |= {processing_date}

        logger.info(f"preparing coarse image for STARS albedo at {cl.place(tile)} on {cl.time(processing_date)}")

        try:
            albedo_coarse_image = generate_albedo_coarse_image(
                date_UTC=processing_date,
                VIIRS_connection=albedo_VIIRS_connection,
                geometry=albedo_coarse_geometry
            )

            albedo_coarse_filename = generate_filename(
                directory=albedo_coarse_directory,
                variable="albedo",
                date_UTC=processing_date,
                tile=tile,
                cell_size=albedo_resolution
            )

            logger.info(f"saving coarse image for STARS albedo at {cl.place(tile)} on {cl.time(processing_date)}: {albedo_coarse_filename}")
            albedo_coarse_image.to_geotiff(albedo_coarse_filename)

            if processing_date >= HLS_start_date:
                logger.info(f"preparing fine image for STARS albedo at {cl.place(tile)} on {cl.time(processing_date)}")

                try:
                    albedo_fine_image = generate_albedo_fine_image(
                        date_UTC=processing_date,
                        tile=tile,
                        HLS_connection=HLS_connection
                    )

                    if calibrate_fine:
                        logger.info(f"calibrating fine image for STARS albedo at {cl.place(tile)} on {cl.time(processing_date)}")
                        albedo_fine_image = calibrate_fine_to_coarse(albedo_fine_image, albedo_coarse_image)

                    albedo_fine_filename = generate_filename(
                        directory=albedo_fine_directory,
                        variable="albedo",
                        date_UTC=processing_date,
                        tile=tile,
                        cell_size=target_resolution
                    )

                    logger.info(f"saving fine image for STARS albedo at {cl.place(tile)} on {cl.time(processing_date)}: {albedo_fine_filename}")
                    albedo_fine_image.to_geotiff(albedo_fine_filename)
                except Exception as e:
                    logger.info(f"HLS is not available on {processing_date}")
        except Exception as e:
            logger.warning(f"unable to produce coarse albedo for date {processing_date}")
            missing_coarse_dates |= {processing_date}

    coarse_latency_dates = [d for d in missing_coarse_dates if (datetime.utcnow().date() - d).days <= VIIRS_GIVEUP_DAYS]

    if len(coarse_latency_dates) > 0:
        raise AncillaryLatency(f"missing coarse dates within {VIIRS_GIVEUP_DAYS}-day window: {', '.join([str(d) for d in coarse_latency_dates])}")

def load_prior(
        tile: str,
        target_resolution: int,
        model_directory: str,
        L2T_STARS_prior_filename: str) -> Prior:
    using_prior = False
    prior_date_UTC = None
    prior_NDVI_filename = None
    prior_NDVI_UQ_filename = None
    prior_NDVI_bias_filename = None
    prior_NDVI_bias_UQ_filename = None
    prior_albedo_filename = None
    prior_albedo_UQ_filename = None
    prior_albedo_bias_filename = None
    prior_albedo_bias_UQ_filename = None

    # check if prior L2T_STARS product is available
    if L2T_STARS_prior_filename is not None and exists(L2T_STARS_prior_filename):
        logger.info(f"loading prior L2T STARS product: {L2T_STARS_prior_filename}")
        L2T_STARS_prior_granule = L2TSTARS(L2T_STARS_prior_filename)
        prior_date_UTC = L2T_STARS_prior_granule.date_UTC
        logger.info(f"prior date: {prior_date_UTC}")
        NDVI_prior_mean = L2T_STARS_prior_granule.NDVI
        NDVI_prior_UQ = L2T_STARS_prior_granule.NDVI_UQ
        albedo_prior_mean = L2T_STARS_prior_granule.albedo
        albedo_prior_UQ = L2T_STARS_prior_granule.albedo_UQ

        prior_tile_date_directory = generate_model_state_tile_date_directory(
            model_directory=model_directory,
            tile=tile,
            date_UTC=prior_date_UTC
        )

        prior_NDVI_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="NDVI",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        NDVI_prior_mean.to_geotiff(prior_NDVI_filename)

        prior_NDVI_UQ_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="NDVI.UQ",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        NDVI_prior_UQ.to_geotiff(prior_NDVI_UQ_filename)

        prior_NDVI_bias_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="NDVI.bias",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        prior_NDVI_bias_UQ_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="NDVI.bias.UQ",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        prior_albedo_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="albedo",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        albedo_prior_mean.to_geotiff(prior_albedo_filename)

        prior_albedo_UQ_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="albedo.UQ",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        albedo_prior_UQ.to_geotiff(prior_albedo_UQ_filename)

        prior_albedo_bias_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="albedo.bias",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        prior_albedo_bias_UQ_filename = generate_filename(
            directory=prior_tile_date_directory,
            variable="albedo.bias.UQ",
            date_UTC=prior_date_UTC,
            tile=tile,
            cell_size=target_resolution
        )

        using_prior = True

    if prior_NDVI_filename is not None and exists(prior_NDVI_filename):
        logger.info(f"prior NDVI ready: {prior_NDVI_filename}")
    else:
        logger.info(f"prior NDVI not found: {prior_NDVI_filename}")
        using_prior = False

    if prior_NDVI_UQ_filename is not None and exists(prior_NDVI_UQ_filename):
        logger.info(f"prior NDVI UQ ready: {prior_NDVI_UQ_filename}")
    else:
        logger.info(f"prior NDVI UQ not found: {prior_NDVI_UQ_filename}")
        using_prior = False

    if prior_NDVI_bias_filename is not None and exists(prior_NDVI_bias_filename):
        logger.info(f"prior NDVI bias ready: {prior_NDVI_bias_filename}")
    else:
        logger.info(f"prior NDVI bias not found: {prior_NDVI_bias_filename}")
        using_prior = False

    if prior_NDVI_bias_UQ_filename is not None and exists(prior_NDVI_bias_UQ_filename):
        logger.info(f"prior NDVI bias UQ ready: {prior_NDVI_bias_UQ_filename}")
    else:
        logger.info(f"prior NDVI bias UQ not found: {prior_NDVI_bias_UQ_filename}")
        using_prior = False

    if prior_albedo_filename is not None and exists(prior_albedo_filename):
        logger.info(f"prior albedo ready: {prior_albedo_filename}")
    else:
        logger.info(f"prior albedo not found: {prior_albedo_filename}")
        using_prior = False

    if prior_albedo_UQ_filename is not None and exists(prior_albedo_UQ_filename):
        logger.info(f"prior albedo UQ ready: {prior_albedo_UQ_filename}")
    else:
        logger.info(f"prior albedo UQ not found: {prior_albedo_UQ_filename}")
        using_prior = False

    if prior_albedo_bias_filename is not None and exists(prior_albedo_bias_filename):
        logger.info(f"prior albedo bias ready: {prior_albedo_bias_filename}")
    else:
        logger.info(f"prior albedo bias not found: {prior_albedo_bias_filename}")
        using_prior = False

    if prior_albedo_bias_UQ_filename is not None and exists(prior_albedo_bias_UQ_filename):
        logger.info(f"prior albedo bias UQ ready: {prior_albedo_bias_UQ_filename}")
    else:
        logger.info(f"prior albedo bias UQ not found: {prior_albedo_bias_UQ_filename}")
        using_prior = False

    prior = Prior(
        using_prior=using_prior,
        prior_date_UTC=prior_date_UTC,
        L2T_STARS_prior_filename=L2T_STARS_prior_filename,
        prior_NDVI_filename=prior_NDVI_filename,
        prior_NDVI_UQ_filename=prior_NDVI_UQ_filename,
        prior_NDVI_bias_filename=prior_NDVI_bias_filename,
        prior_NDVI_bias_UQ_filename=prior_NDVI_bias_UQ_filename,
        prior_albedo_filename=prior_albedo_filename,
        prior_albedo_UQ_filename=prior_albedo_UQ_filename,
        prior_albedo_bias_filename=prior_albedo_bias_filename,
        prior_albedo_bias_UQ_filename=prior_albedo_bias_UQ_filename
    )

    return prior

def process_STARS_product(
        tile: str,
        date_UTC: date,
        time_UTC: datetime,
        build: str,
        product_counter: int,
        HLS_start_date: date,
        HLS_end_date: date,
        VIIRS_start_date: date,
        VIIRS_end_date: date,
        NDVI_resolution: int,
        albedo_resolution: int,
        target_resolution: int,
        working_directory: str,
        model_directory: str,
        input_staging_directory: str,
        L2T_STARS_granule_directory: str,
        L2T_STARS_zip_filename: str,
        L2T_STARS_browse_filename: str,
        metadata: dict,
        prior: Prior,
        HLS_connection: HLS2CMR,
        NDVI_VIIRS_connection: VIIRSDownloaderNDVI,
        albedo_VIIRS_connection: VIIRSDownloaderAlbedo,
        using_prior: bool = False,
        calibrate_fine: bool = True,
        remove_input_staging: bool = True,
        remove_prior: bool = True,
        remove_posterior: bool = True):
    NDVI_coarse_geometry = HLS_connection.grid(tile=tile, cell_size=NDVI_resolution)
    albedo_coarse_geometry = HLS_connection.grid(tile=tile, cell_size=albedo_resolution)

    logger.info(f"processing the L2T_STARS product at tile {tile} for date {date_UTC}")

    NDVI_coarse_directory = generate_NDVI_coarse_directory(
        input_staging_directory=input_staging_directory,
        tile=tile
    )

    logger.info(f"staging coarse NDVI images: {NDVI_coarse_directory}")

    NDVI_fine_directory = generate_NDVI_fine_directory(
        input_staging_directory=input_staging_directory,
        tile=tile
    )

    logger.info(f"staging fine NDVI images: {NDVI_fine_directory}")

    albedo_coarse_directory = generate_albedo_coarse_directory(
        input_staging_directory=input_staging_directory,
        tile=tile
    )

    logger.info(f"staging coarse albedo images: {albedo_coarse_directory}")

    albedo_fine_directory = generate_albedo_fine_directory(
        input_staging_directory=input_staging_directory,
        tile=tile
    )

    logger.info(f"staging fine albedo images: {albedo_fine_directory}")

    posterior_tile_date_directory = generate_model_state_tile_date_directory(
        model_directory=model_directory,
        tile=tile,
        date_UTC=date_UTC
    )

    logger.info(f"posterior directory: {posterior_tile_date_directory}")

    generate_STARS_inputs(
        tile=tile,
        date_UTC=date_UTC,
        HLS_start_date=HLS_start_date,
        HLS_end_date=HLS_end_date,
        VIIRS_start_date=VIIRS_start_date,
        VIIRS_end_date=VIIRS_end_date,
        NDVI_resolution=NDVI_resolution,
        albedo_resolution=albedo_resolution,
        target_resolution=target_resolution,
        NDVI_coarse_geometry=NDVI_coarse_geometry,
        albedo_coarse_geometry=albedo_coarse_geometry,
        working_directory=working_directory,
        NDVI_coarse_directory=NDVI_coarse_directory,
        NDVI_fine_directory=NDVI_fine_directory,
        albedo_coarse_directory=albedo_coarse_directory,
        albedo_fine_directory=albedo_fine_directory,
        HLS_connection=HLS_connection,
        NDVI_VIIRS_connection=NDVI_VIIRS_connection,
        albedo_VIIRS_connection=albedo_VIIRS_connection,
        calibrate_fine=calibrate_fine
    )

    posterior_NDVI_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="NDVI",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    logger.info(f"posterior NDVI file: {posterior_NDVI_filename}")

    posterior_NDVI_UQ_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="NDVI.UQ",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    logger.info(f"posterior NDVI UQ file: {posterior_NDVI_UQ_filename}")

    posterior_NDVI_bias_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="NDVI.bias",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    logger.info(f"posterior NDVI bias file: {posterior_NDVI_bias_filename}")

    posterior_NDVI_bias_UQ_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="NDVI.bias.UQ",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    logger.info(f"posterior NDVI bias UQ file: {posterior_NDVI_bias_UQ_filename}")

    if using_prior:
        process_julia_data_fusion(
            tile=tile,
            coarse_cell_size=NDVI_resolution,
            fine_cell_size=target_resolution,
            VIIRS_start_date=VIIRS_start_date,
            VIIRS_end_date=VIIRS_end_date,
            HLS_start_date=HLS_start_date,
            HLS_end_date=HLS_end_date,
            coarse_directory=NDVI_coarse_directory,
            fine_directory=NDVI_fine_directory,
            posterior_filename=posterior_NDVI_filename,
            posterior_UQ_filename=posterior_NDVI_UQ_filename,
            posterior_bias_filename=posterior_NDVI_bias_filename,
            posterior_bias_UQ_filename=posterior_NDVI_bias_UQ_filename,
            prior_filename=prior.prior_NDVI_filename,
            prior_UQ_filename=prior.prior_NDVI_UQ_filename,
            prior_bias_filename=prior.prior_NDVI_bias_filename,
            prior_bias_UQ_filename=prior.prior_NDVI_bias_UQ_filename
        )
    else:
        process_julia_data_fusion(
            tile=tile,
            coarse_cell_size=NDVI_resolution,
            fine_cell_size=target_resolution,
            VIIRS_start_date=VIIRS_start_date,
            VIIRS_end_date=VIIRS_end_date,
            HLS_start_date=HLS_start_date,
            HLS_end_date=HLS_end_date,
            coarse_directory=NDVI_coarse_directory,
            fine_directory=NDVI_fine_directory,
            posterior_filename=posterior_NDVI_filename,
            posterior_UQ_filename=posterior_NDVI_UQ_filename,
            posterior_bias_filename=posterior_NDVI_bias_filename,
            posterior_bias_UQ_filename=posterior_NDVI_bias_UQ_filename
        )

    NDVI = Raster.open(posterior_NDVI_filename)
    NDVI_UQ = Raster.open(posterior_NDVI_UQ_filename)

    posterior_albedo_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="albedo",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    posterior_albedo_UQ_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="albedo.UQ",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    posterior_albedo_bias_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="albedo.bias",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    posterior_albedo_bias_UQ_filename = generate_filename(
        directory=posterior_tile_date_directory,
        variable="albedo.bias.UQ",
        date_UTC=date_UTC,
        tile=tile,
        cell_size=target_resolution
    )

    if using_prior:
        process_julia_data_fusion(
            tile=tile,
            coarse_cell_size=albedo_resolution,
            fine_cell_size=target_resolution,
            VIIRS_start_date=VIIRS_start_date,
            VIIRS_end_date=VIIRS_end_date,
            HLS_start_date=HLS_start_date,
            HLS_end_date=HLS_end_date,
            coarse_directory=albedo_coarse_directory,
            fine_directory=albedo_fine_directory,
            posterior_filename=posterior_albedo_filename,
            posterior_UQ_filename=posterior_albedo_UQ_filename,
            posterior_bias_filename=posterior_albedo_bias_filename,
            posterior_bias_UQ_filename=posterior_albedo_bias_UQ_filename,
            prior_filename=prior.prior_albedo_filename,
            prior_UQ_filename=prior.prior_albedo_UQ_filename,
            prior_bias_filename=prior.prior_albedo_bias_filename,
            prior_bias_UQ_filename=prior.prior_albedo_bias_UQ_filename,
        )
    else:
        process_julia_data_fusion(
            tile=tile,
            coarse_cell_size=albedo_resolution,
            fine_cell_size=target_resolution,
            VIIRS_start_date=VIIRS_start_date,
            VIIRS_end_date=VIIRS_end_date,
            HLS_start_date=HLS_start_date,
            HLS_end_date=HLS_end_date,
            coarse_directory=albedo_coarse_directory,
            fine_directory=albedo_fine_directory,
            posterior_filename=posterior_albedo_filename,
            posterior_UQ_filename=posterior_albedo_UQ_filename,
            posterior_bias_filename=posterior_albedo_bias_filename,
            posterior_bias_UQ_filename=posterior_albedo_bias_UQ_filename
        )

    albedo = Raster.open(posterior_albedo_filename)
    albedo_UQ = Raster.open(posterior_albedo_UQ_filename)

    if NDVI is None:
        raise BlankOutput("unable to generate STARS NDVI")

    if NDVI_UQ is None:
        raise BlankOutput("unable to generate STARS NDVI UQ")

    if albedo is None:
        raise BlankOutput("unable to generate STARS albedo")

    if albedo_UQ is None:
        raise BlankOutput("unable to generate STARS albedo UQ")

    granule = L2TSTARS(
        product_location=L2T_STARS_granule_directory,
        tile=tile,
        time_UTC=time_UTC,
        build=build,
        process_count=product_counter
    )

    granule.add_layer("NDVI", NDVI, cmap=NDVI_COLORMAP)
    granule.add_layer("NDVI-UQ", NDVI_UQ, cmap="jet")
    granule.add_layer("albedo", albedo, cmap=ALBEDO_COLORMAP)
    granule.add_layer("albedo-UQ", albedo_UQ, cmap="jet")

    metadata["StandardMetadata"]["LocalGranuleID"] = basename(
        L2T_STARS_zip_filename)
    metadata["StandardMetadata"]["SISName"] = "Level 2 STARS Product Specification Document"
    granule.write_metadata(metadata)

    logger.info(
        f"writing L2T STARS product zip: {cl.file(L2T_STARS_zip_filename)}")
    granule.write_zip(L2T_STARS_zip_filename)
    logger.info(
        f"writing L2T STARS browse image: {cl.file(L2T_STARS_browse_filename)}")
    granule.write_browse_image(
        PNG_filename=L2T_STARS_browse_filename)
    logger.info(
        f"removing L2T STARS tile granule directory: {cl.dir(L2T_STARS_granule_directory)}")
    shutil.rmtree(L2T_STARS_granule_directory)

    if exists(L2T_STARS_browse_filename):
        logger.info(
            f"found L2T STARS browse image: {cl.file(L2T_STARS_zip_filename)}")
    else:
        logger.info(
            f"generating L2T STARS browse image: {cl.file(L2T_STARS_browse_filename)}")
        granule = L2TSTARS(L2T_STARS_zip_filename)
        granule.write_browse_image(
            PNG_filename=L2T_STARS_browse_filename)

    logger.info(f"re-writing posterior NDVI: {posterior_NDVI_filename}")
    Raster.open(posterior_NDVI_filename, cmap=NDVI_COLORMAP).to_geotiff(posterior_NDVI_filename)
    logger.info(f"re-writing posterior NDVI UQ: {posterior_NDVI_UQ_filename}")
    Raster.open(posterior_NDVI_UQ_filename, cmap=NDVI_COLORMAP).to_geotiff(posterior_NDVI_UQ_filename)
    logger.info(f"re-writing posterior NDVI bias: {posterior_NDVI_bias_filename}")
    Raster.open(posterior_NDVI_bias_filename, cmap=NDVI_COLORMAP).to_geotiff(posterior_NDVI_bias_filename)
    logger.info(f"re-writing posterior NDVI bias UQ: {posterior_NDVI_bias_UQ_filename}")
    Raster.open(posterior_NDVI_bias_UQ_filename, cmap=NDVI_COLORMAP).to_geotiff(posterior_NDVI_bias_UQ_filename)

    logger.info(f"re-writing posterior albedo: {posterior_albedo_filename}")
    Raster.open(posterior_albedo_filename, cmap=ALBEDO_COLORMAP).to_geotiff(posterior_albedo_filename)
    logger.info(f"re-writing posterior albedo UQ: {posterior_albedo_UQ_filename}")
    Raster.open(posterior_albedo_UQ_filename, cmap=ALBEDO_COLORMAP).to_geotiff(posterior_albedo_UQ_filename)
    logger.info(f"re-writing posterior albedo bias: {posterior_albedo_bias_filename}")
    Raster.open(posterior_albedo_bias_filename, cmap=ALBEDO_COLORMAP).to_geotiff(posterior_albedo_bias_filename)
    logger.info(f"re-writing posterior albedo bias UQ: {posterior_albedo_bias_UQ_filename}")
    Raster.open(posterior_albedo_bias_UQ_filename, cmap=ALBEDO_COLORMAP).to_geotiff(posterior_albedo_bias_UQ_filename)

    if remove_input_staging:
        logger.info(f"removing input staging directory: {input_staging_directory}")
        shutil.rmtree(input_staging_directory)

    if using_prior and remove_prior:
        if exists(prior.prior_NDVI_filename):
            logger.info(f"removing NDVI prior: {prior.prior_NDVI_filename}")
            remove(prior.prior_NDVI_filename)

        if exists(prior.prior_NDVI_UQ_filename):
            logger.info(f"removing NDVI UQ prior: {prior.prior_NDVI_UQ_filename}")
            remove(prior.prior_NDVI_UQ_filename)

        if exists(prior.prior_NDVI_bias_filename):
            logger.info(f"removing NDVI bias prior: {prior.prior_NDVI_bias_filename}")
            remove(prior.prior_NDVI_bias_filename)

        if exists(prior.prior_NDVI_bias_UQ_filename):
            logger.info(f"removing NDVI bias UQ prior: {prior.prior_NDVI_bias_UQ_filename}")
            remove(prior.prior_NDVI_bias_UQ_filename)

        if exists(prior.prior_albedo_filename):
            logger.info(f"removing albedo prior: {prior.prior_albedo_filename}")
            remove(prior.prior_albedo_filename)

        if exists(prior.prior_albedo_UQ_filename):
            logger.info(f"removing albedo UQ prior: {prior.prior_albedo_UQ_filename}")
            remove(prior.prior_albedo_UQ_filename)

        if exists(prior.prior_albedo_bias_filename):
            logger.info(f"removing albedo bias prior: {prior.prior_albedo_bias_filename}")
            remove(prior.prior_albedo_bias_filename)

        if exists(prior.prior_albedo_bias_UQ_filename):
            logger.info(f"removing albedo bias UQ prior: {prior.prior_albedo_bias_UQ_filename}")
            remove(prior.prior_albedo_bias_UQ_filename)

    if remove_posterior:
        if exists(posterior_NDVI_filename):
            logger.info(f"removing NDVI posterior: {posterior_NDVI_filename}")
            remove(posterior_NDVI_filename)

        if exists(posterior_NDVI_UQ_filename):
            logger.info(f"removing NDVI UQ posterior: {posterior_NDVI_UQ_filename}")
            remove(posterior_NDVI_UQ_filename)

        if exists(posterior_albedo_filename):
            logger.info(f"removing albedo posterior: {posterior_albedo_filename}")
            remove(posterior_albedo_filename)

        if exists(posterior_albedo_UQ_filename):
            logger.info(f"removing albedo UQ posterior: {posterior_albedo_UQ_filename}")
            remove(posterior_albedo_UQ_filename)

def L2T_STARS(
        runconfig_filename: str,
        date_UTC: Union[date, str] = None,
        spinup_days: int = DEFAULT_SPINUP_DAYS,
        target_resolution: int = DEFAULT_TARGET_RESOLUTION,
        NDVI_resolution: int = DEFAULT_NDVI_RESOLUTION,
        albedo_resolution: int = DEFAULT_ALBEDO_RESOLUTION,
        use_VNP43NRT: bool = DEFAULT_USE_VNP43NRT,
        calibrate_fine: bool = DEFAULT_CALIBRATE_FINE,
        sources_only: bool = False,
        remove_input_staging: bool = True,
        remove_prior: bool = True,
        remove_posterior: bool = True) -> int:
    """
    ECOSTRESS Collection 2 L2G L2T LSTE PGE
    :param runconfig_filename: filename for XML run-config
    :param log_filename: filename for logger output
    :return: exit code number
    """
    exit_code = SUCCESS_EXIT_CODE

    try:
        runconfig = L2TSTARSConfig(runconfig_filename)
        working_directory = runconfig.working_directory
        granule_ID = runconfig.granule_ID
        log_filename = join(working_directory, "log", f"{granule_ID}.log")
        cl.configure(filename=log_filename)

        logger.info(f"L2T_STARS PGE ({cl.val(ECOSTRESS.PGEVersion)})")
        logger.info(f"L2T_STARS run-config: {cl.file(runconfig_filename)}")

        logger.info(f"granule ID: {granule_ID}")

        L2T_STARS_granule_directory = runconfig.L2T_STARS_granule_directory
        logger.info(f"granule directory: {cl.dir(L2T_STARS_granule_directory)}")
        L2T_STARS_zip_filename = runconfig.L2T_STARS_zip_filename
        logger.info(f"zip filename: {cl.file(L2T_STARS_zip_filename)}")
        L2T_STARS_browse_filename = runconfig.L2T_STARS_browse_filename
        logger.info(f"browse filename: " + cl.file(L2T_STARS_browse_filename))

        if exists(L2T_STARS_zip_filename) and exists(L2T_STARS_browse_filename):
            logger.info(f"found L2T STARS file: {L2T_STARS_zip_filename}")
            logger.info(
                f"found L2T STARS preview: {L2T_STARS_browse_filename}")
            return SUCCESS_EXIT_CODE

        logger.info(f"working_directory: {cl.dir(working_directory)}")
        logger.info(f"log: {cl.file(log_filename)}")

        input_staging_directory = join(working_directory, "input_staging")

        sources_directory = runconfig.sources_directory
        logger.info(f"source directory: {cl.dir(sources_directory)}")
        indices_directory = runconfig.indices_directory
        logger.info(f"indices directory: {cl.dir(indices_directory)}")
        model_directory = runconfig.model_directory
        logger.info(f"model directory: {cl.dir(model_directory)}")
        output_directory = runconfig.output_directory
        logger.info(f"output directory: {cl.dir(output_directory)}")
        tile = runconfig.tile
        logger.info(f"tile: {cl.val(tile)}")
        build = runconfig.build
        logger.info(f"build: {cl.val(build)}")
        product_counter = runconfig.product_counter
        logger.info(f"product counter: {cl.val(product_counter)}")
        L2T_LSTE_filename = runconfig.L2T_LSTE_filename
        logger.info(f"L2T LSTE file: {cl.file(L2T_LSTE_filename)}")

        if not exists(L2T_LSTE_filename):
            raise InputFilesInaccessible(
                f"L2T LSTE file does not exist: {L2T_LSTE_filename}")

        L2T_granule = L2TLSTE(L2T_LSTE_filename)
        geometry = L2T_granule.geometry
        metadata = L2T_granule.metadata_dict
        metadata["StandardMetadata"]["PGEName"] = "L2T_STARS"

        short_name = L2T_STARS_SHORT_NAME
        logger.info(f"L2T STARS short name: {cl.val(short_name)}")
        metadata["StandardMetadata"]["ShortName"] = short_name

        long_name = L2T_STARS_LONG_NAME
        logger.info(f"L2T STARS long name: {cl.val(long_name)}")
        metadata["StandardMetadata"]["LongName"] = long_name

        metadata["StandardMetadata"]["AncillaryInputPointer"] = "HLS,VIIRS"
        del (metadata["ProductMetadata"]["AncillaryNWP"])
        del (metadata["ProductMetadata"]["NWPSource"])

        time_UTC = L2T_granule.time_UTC
        logger.info(f"ECOSTRESS overpass time: {cl.time(f'{time_UTC:%Y-%m-%d %H:%M:%S} UTC')}")

        if date_UTC is None:
            date_UTC = L2T_granule.date_UTC
            logger.info(f"ECOSTRESS overpass date: {cl.time(f'{date_UTC:%Y-%m-%d} UTC')}")
        else:
            logger.warning(f"over-riding target date: {date_UTC}")

            if isinstance(date_UTC, str):
                date_UTC = parser.parse(date_UTC).date()

        L2T_STARS_prior_filename = runconfig.L2T_STARS_prior_filename

        prior = load_prior(
            tile=tile,
            target_resolution=target_resolution,
            model_directory=model_directory,
            L2T_STARS_prior_filename=L2T_STARS_prior_filename
        )

        using_prior = prior.using_prior
        prior_date_UTC = prior.prior_date_UTC

        products_directory = join(working_directory, DEFAULT_STARS_PRODUCTS_DIRECTORY)
        logger.info(f"STARS products directory: {cl.dir(products_directory)}")
        HLS_download_directory = join(sources_directory, DEFAULT_HLS_DOWNLOAD_DIRECTORY)
        logger.info(f"HLS download directory: {cl.dir(HLS_download_directory)}")
        HLS_products_directory = join(sources_directory, DEFAULT_HLS_PRODUCTS_DIRECTORY)
        logger.info(f"HLS products directory: {cl.dir(HLS_products_directory)}")
        VIIRS_download_directory = join(sources_directory, DEFAULT_VIIRS_DOWNLOAD_DIRECTORY)
        logger.info(f"VIIRS download directory: {cl.dir(VIIRS_download_directory)}")
        VIIRS_products_directory = join(sources_directory, DEFAULT_VIIRS_PRODUCTS_DIRECTORY)
        logger.info(f"VIIRS products directory: {cl.dir(VIIRS_products_directory)}")
        VIIRS_mosaic_directory = join(sources_directory, DEFAUL_VIIRS_MOSAIC_DIRECTORY)
        logger.info(f"VIIRS mosaic directory: {cl.dir(VIIRS_mosaic_directory)}")
        GEOS5FP_download_directory = join(sources_directory, DEFAULT_GEOS5FP_DOWNLOAD_DIRECTORY)
        logger.info(f"GEOS-5 FP download directory: {cl.dir(GEOS5FP_download_directory)}")
        GEOS5FP_products_directory = join(sources_directory, DEFAULT_GEOS5FP_PRODUCTS_DIRECTORY)
        logger.info(f"GEOS-5 FP products directory: {cl.dir(GEOS5FP_products_directory)}")
        VNP09GA_products_directory = join(sources_directory, DEFAULT_VNP09GA_PRODUCTS_DIRECTORY)
        logger.info(f"VNP09GA products directory: {cl.dir(VNP09GA_products_directory)}")
        VNP43NRT_products_directory = join(sources_directory, DEFAULT_VNP43NRT_PRODUCTS_DIRECTORY)
        logger.info(f"VNP43NRT products directory: {cl.dir(VNP43NRT_products_directory)}")

        if exists(L2T_STARS_zip_filename):
            logger.info(
                f"found L2T STARS product zip: {cl.file(L2T_STARS_zip_filename)}")
            return exit_code

        logger.info(f"connecting to CMR Search server: {CMR_SEARCH_URL}")
        try:
            HLS_connection = HLS2CMR(
                working_directory=working_directory,
                download_directory=HLS_download_directory,
                products_directory=HLS_products_directory,
                target_resolution=target_resolution
            )
        except CMRServerUnreachable as e:
            logger.exception(e)
            raise AncillaryServerUnreachable(
                f"unable to connect to CMR Search server: {CMR_SEARCH_URL}")

        if not HLS_connection.tile_grid.land(tile=tile):
            raise LandFilter(f"Sentinel tile {tile} is not on land")

        if use_VNP43NRT:
            try:
                NDVI_VIIRS_connection = VNP43NRT(
                    working_directory=working_directory,
                    download_directory=VIIRS_download_directory,
                    mosaic_directory=VIIRS_mosaic_directory,
                    GEOS5FP_download=GEOS5FP_download_directory,
                    GEOS5FP_products=GEOS5FP_products_directory,
                    VNP09GA_directory=VNP09GA_products_directory,
                    VNP43NRT_directory=VNP43NRT_products_directory
                )

                albedo_VIIRS_connection = VNP43NRT(
                    working_directory=working_directory,
                    download_directory=VIIRS_download_directory,
                    mosaic_directory=VIIRS_mosaic_directory,
                    GEOS5FP_download=GEOS5FP_download_directory,
                    GEOS5FP_products=GEOS5FP_products_directory,
                    VNP09GA_directory=VNP09GA_products_directory,
                    VNP43NRT_directory=VNP43NRT_products_directory
                )
            except CMRServerUnreachable as e:
                logger.exception(e)
                raise AncillaryServerUnreachable(
                    f"unable to connect to CMR search server")
        else:
            try:
                NDVI_VIIRS_connection = VNP43IA4(
                    working_directory=working_directory,
                    download_directory=VIIRS_download_directory,
                    products_directory=VIIRS_products_directory,
                    mosaic_directory=VIIRS_mosaic_directory
                )

                albedo_VIIRS_connection = VNP43MA3(
                    working_directory=working_directory,
                    download_directory=VIIRS_download_directory,
                    products_directory=VIIRS_products_directory,
                    mosaic_directory=VIIRS_mosaic_directory
                )
            except LPDAACServerUnreachable as e:
                logger.exception(e)
                raise AncillaryServerUnreachable(
                    f"unable to connect to VIIRS server")

        end_date = date_UTC
        # the start date of the BRDF-corrected VIIRS coarse time-series is 10 days before the target date
        VIIRS_start_date = end_date - timedelta(days=spinup_days)
        # to produce that first BRDF-corrected image, we need VNP09GA starting 16 days prior to the first coarse date
        VIIRS_download_start_date = VIIRS_start_date - timedelta(days=16)
        VIIRS_end_date = end_date

        # define start date of HLS fine image input time-series
        if using_prior and prior_date_UTC >= VIIRS_start_date:
            # if we're using a prior, HLS inputs begin the day after the prior
            HLS_start_date = prior_date_UTC + timedelta(days=1)
        else:
            # if we're initializing, HLS inputs begin on the same day as the VIIRS inputs
            HLS_start_date = VIIRS_start_date

        HLS_end_date = end_date

        logger.info(
            f"processing STARS HLS-VIIRS NDVI and albedo for tile {cl.place(tile)} from {cl.time(VIIRS_start_date)} to {cl.time(end_date)}")

        try:
            HLS_listing = HLS_connection.listing(
                tile=tile,
                start_UTC=HLS_start_date,
                end_UTC=HLS_end_date
            )

        except HLSTileNotAvailable as e:
            logger.exception(e)

            raise LandFilter(
                f"Sentinel tile {tile} cannot be processed")
        except Exception as e:
            logger.exception(e)
            raise AncillaryServerUnreachable(
                f"unable to scan Harmonized Landsat Sentinel server: {HLS_connection.remote}")

        missing_sentinel_dates = HLS_listing[HLS_listing.sentinel == "missing"].date_UTC

        if len(missing_sentinel_dates) > 0:
            raise AncillaryLatency(
                f"HLS Sentinel is not yet available at tile {tile} for dates: {', '.join(missing_sentinel_dates)}")

        sentinel_listing = HLS_listing[~pd.isna(HLS_listing.sentinel)][["date_UTC", "sentinel"]]

        logger.info(f"HLS Sentinel is available on {cl.val(len(sentinel_listing))} dates:")
        for i, (date_UTC, sentinel_granule) in sentinel_listing.iterrows():
            sentinel_filename = sentinel_granule["meta"]["native-id"]
            logger.info(f"* {cl.time(date_UTC)}: {cl.file(sentinel_filename)}")

        missing_landsat_dates = HLS_listing[HLS_listing.landsat == "missing"].date_UTC

        if len(missing_landsat_dates) > 0:
            raise AncillaryLatency(
                f"HLS landsat is not yet available at tile {tile} for dates: {', '.join(missing_landsat_dates)}")

        landsat_listing = HLS_listing[~pd.isna(HLS_listing.landsat)][["date_UTC", "landsat"]]

        logger.info(f"HLS landsat is available on {cl.val(len(landsat_listing))} dates:")
        for i, (date_UTC, landsat_granule) in landsat_listing.iterrows():
            landsat_filename = landsat_granule["meta"]["native-id"]
            logger.info(f"* {cl.time(date_UTC)}: {cl.file(landsat_filename)}")

        if sources_only:
            retrieve_STARS_sources(
                tile=tile,
                geometry=geometry,
                HLS_start_date=HLS_start_date,
                HLS_end_date=HLS_end_date,
                VIIRS_start_date=VIIRS_download_start_date,
                VIIRS_end_date=VIIRS_end_date,
                HLS_connection=HLS_connection,
                VIIRS_connection=NDVI_VIIRS_connection,
            )

            NDVI_coarse_geometry = HLS_connection.grid(tile=tile, cell_size=NDVI_resolution)
            albedo_coarse_geometry = HLS_connection.grid(tile=tile, cell_size=albedo_resolution)

            logger.info(f"processing the L2T_STARS product at tile {tile} for date {date_UTC}")

            NDVI_coarse_directory = generate_NDVI_coarse_directory(
                input_staging_directory=input_staging_directory,
                tile=tile
            )

            logger.info(f"staging coarse NDVI images: {NDVI_coarse_directory}")

            NDVI_fine_directory = generate_NDVI_fine_directory(
                input_staging_directory=input_staging_directory,
                tile=tile
            )

            logger.info(f"staging fine NDVI images: {NDVI_fine_directory}")

            albedo_coarse_directory = generate_albedo_coarse_directory(
                input_staging_directory=input_staging_directory,
                tile=tile
            )

            logger.info(f"staging coarse albedo images: {albedo_coarse_directory}")

            albedo_fine_directory = generate_albedo_fine_directory(
                input_staging_directory=input_staging_directory,
                tile=tile
            )

            logger.info(f"staging fine albedo images: {albedo_fine_directory}")

            posterior_tile_date_directory = generate_model_state_tile_date_directory(
                model_directory=model_directory,
                tile=tile,
                date_UTC=date_UTC
            )

            logger.info(f"posterior directory: {posterior_tile_date_directory}")

            generate_STARS_inputs(
                tile=tile,
                date_UTC=date_UTC,
                HLS_start_date=HLS_start_date,
                HLS_end_date=HLS_end_date,
                VIIRS_start_date=VIIRS_start_date,
                VIIRS_end_date=VIIRS_end_date,
                NDVI_resolution=NDVI_resolution,
                albedo_resolution=albedo_resolution,
                target_resolution=target_resolution,
                NDVI_coarse_geometry=NDVI_coarse_geometry,
                albedo_coarse_geometry=albedo_coarse_geometry,
                working_directory=working_directory,
                NDVI_coarse_directory=NDVI_coarse_directory,
                NDVI_fine_directory=NDVI_fine_directory,
                albedo_coarse_directory=albedo_coarse_directory,
                albedo_fine_directory=albedo_fine_directory,
                HLS_connection=HLS_connection,
                NDVI_VIIRS_connection=NDVI_VIIRS_connection,
                albedo_VIIRS_connection=albedo_VIIRS_connection,
                calibrate_fine=calibrate_fine
            )
        else:
            process_STARS_product(
                tile=tile,
                date_UTC=date_UTC,
                time_UTC=time_UTC,
                build=build,
                product_counter=product_counter,
                HLS_start_date=HLS_start_date,
                HLS_end_date=HLS_end_date,
                VIIRS_start_date=VIIRS_start_date,
                VIIRS_end_date=VIIRS_end_date,
                NDVI_resolution=NDVI_resolution,
                albedo_resolution=albedo_resolution,
                target_resolution=target_resolution,
                working_directory=working_directory,
                model_directory=model_directory,
                input_staging_directory=input_staging_directory,
                L2T_STARS_granule_directory=L2T_STARS_granule_directory,
                L2T_STARS_zip_filename=L2T_STARS_zip_filename,
                L2T_STARS_browse_filename=L2T_STARS_browse_filename,
                metadata=metadata,
                prior=prior,
                HLS_connection=HLS_connection,
                NDVI_VIIRS_connection=NDVI_VIIRS_connection,
                albedo_VIIRS_connection=albedo_VIIRS_connection,
                using_prior=using_prior,
                calibrate_fine=calibrate_fine,
                remove_input_staging=remove_input_staging,
                remove_prior=remove_prior,
                remove_posterior=remove_posterior
            )

    except (ConnectionError, urllib.error.HTTPError, CMRServerUnreachable) as exception:
        logger.exception(exception)
        exit_code = ANCILLARY_SERVER_UNREACHABLE

    except DownloadFailed as exception:
        logger.exception(exception)
        exit_code = DOWNLOAD_FAILED

    except HLSBandNotAcquired as exception:
        logger.exception(exception)
        exit_code = DOWNLOAD_FAILED

    except HLSNotAvailable as exception:
        logger.exception(exception)
        exit_code = LAND_FILTER

    except (HLSSentinelMissing, HLSLandsatMissing) as exception:
        logger.exception(exception)
        exit_code = ANCILLARY_LATENCY

    except ECOSTRESSExitCodeException as exception:
        logger.exception(exception)
        exit_code = exception.exit_code

    except Exception as exception:
        logger.exception(exception)
        exit_code = UNCLASSIFIED_FAILURE_EXIT_CODE

    logger.info(f"L2T_STARS exit code: {exit_code}")

    return exit_code


def main(argv=sys.argv):
    if len(argv) == 1 or "--version" in argv:
        print(f"L2T_STARS PGE ({ECOSTRESS.PGEVersion})")
        print(f"usage: L2T_STARS RunConfig.xml")

        if "--version" in argv:
            return SUCCESS_EXIT_CODE
        else:
            return RUNCONFIG_FILENAME_NOT_SUPPLIED

    runconfig_filename = str(argv[1])
    exit_code = L2T_STARS(runconfig_filename=runconfig_filename)

    return exit_code


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
