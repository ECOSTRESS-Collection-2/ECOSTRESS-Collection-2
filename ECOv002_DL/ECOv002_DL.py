import logging
import socket
import sys
from datetime import datetime
from os import makedirs
from os.path import join, abspath, dirname, expanduser, exists, splitext, basename
from shutil import which
from typing import List
from uuid import uuid4

import ECOSTRESS
import colored_logging as cl
from ECOSTRESS.L2_LSTE import L2GLSTE, L2TLSTE
from ECOSTRESS.exit_codes import LAND_FILTER, SUCCESS_EXIT_CODE, RUNCONFIG_FILENAME_NOT_SUPPLIED, MissingRunConfigValue, \
    ECOSTRESSExitCodeException, UnableToParseRunConfig
from ECOSTRESS.runconfig import ECOSTRESSRunConfig, read_runconfig
from L2T_STARS import DEFAULT_STARS_SOURCES_DIRECTORY, L2T_STARS, generate_L2T_STARS_runconfig, \
    DEFAULT_STARS_INDICES_DIRECTORY, DEFAULT_STARS_MODEL_DIRECTORY
from sentinel_tile_grid import SentinelTileGrid

with open(join(abspath(dirname(__file__)), "version.txt")) as f:
    version = f.read()

__version__ = version

logger = logging.getLogger(__name__)

ECOv002_DL_TEMPLATE = join(abspath(dirname(__file__)), "ECOv002_DL.xml")
DEFAULT_BUILD = "0700"


def generate_downloader_runconfig(
        L2G_LSTE_filename: str,
        L2T_LSTE_filenames: List[str],
        orbit: int = None,
        scene: int = None,
        working_directory: str = None,
        L2T_STARS_sources_directory: str = None,
        L2T_STARS_indices_directory: str = None,
        L2T_STARS_model_directory: str = None,
        executable_filename: str = None,
        runconfig_filename: str = None,
        log_filename: str = None,
        build: str = None,
        processing_node: str = None,
        production_datetime: datetime = None,
        job_ID: str = None,
        instance_ID: str = None,
        product_counter: int = None,
        template_filename: str = None) -> str:
    L2G_LSTE_filename = abspath(expanduser(L2G_LSTE_filename))

    if not exists(L2G_LSTE_filename):
        raise IOError(f"L2G LSTE file not found: {L2G_LSTE_filename}")

    logger.info(f"L2G LSTE file: {cl.file(L2G_LSTE_filename)}")
    source_granule_ID = splitext(basename(L2G_LSTE_filename))[0]
    logger.info(f"source granule ID: {cl.name(source_granule_ID)}")

    if orbit is None:
        orbit = int(source_granule_ID.split("_")[-5])

    logger.info(f"orbit: {cl.val(orbit)}")

    if scene is None:
        scene = int(source_granule_ID.split("_")[-4])

    logger.info(f"scene: {cl.val(scene)}")

    if template_filename is None:
        template_filename = ECOv002_DL_TEMPLATE

    template_filename = abspath(expanduser(template_filename))

    run_ID = f"ECOv002_DL_{orbit:05d}_{scene:05d}"

    if runconfig_filename is None:
        runconfig_filename = join(working_directory, "runconfig", f"{run_ID}.xml")

    runconfig_filename = abspath(expanduser(runconfig_filename))

    if working_directory is None:
        working_directory = run_ID

    working_directory = abspath(expanduser(working_directory))

    if L2T_STARS_sources_directory is None:
        L2T_STARS_sources_directory = join(working_directory, DEFAULT_STARS_SOURCES_DIRECTORY)

    L2T_STARS_sources_directory = abspath(expanduser(L2T_STARS_sources_directory))

    if L2T_STARS_indices_directory is None:
        L2T_STARS_indices_directory = join(working_directory, DEFAULT_STARS_INDICES_DIRECTORY)

    L2T_STARS_indices_directory = abspath(expanduser(L2T_STARS_indices_directory))

    if L2T_STARS_model_directory is None:
        L2T_STARS_model_directory = join(working_directory, DEFAULT_STARS_MODEL_DIRECTORY)

    L2T_STARS_model_directory = abspath(expanduser(L2T_STARS_model_directory))

    if executable_filename is None:
        executable_filename = which("ECOv002_DL")

    if executable_filename is None:
        executable_filename = "ECOv002_DL"

    if log_filename is None:
        log_filename = join(working_directory, f"{run_ID}.log")

    log_filename = abspath(expanduser(log_filename))

    if build is None:
        build = DEFAULT_BUILD

    if processing_node is None:
        processing_node = socket.gethostname()

    if production_datetime is None:
        production_datetime = datetime.utcnow()

    if isinstance(production_datetime, datetime):
        production_datetime = str(production_datetime)

    if job_ID is None:
        job_ID = production_datetime

    if instance_ID is None:
        instance_ID = str(uuid4())

    if product_counter is None:
        product_counter = 1

    working_directory = abspath(expanduser(working_directory))

    logger.info(f"generating run-config for orbit {cl.val(orbit)} scene {cl.val(scene)}")
    logger.info(f"loading ECOv002_DL template: {cl.file(template_filename)}")

    with open(template_filename, "r") as file:
        template = file.read()

    logger.info(f"orbit: {cl.val(orbit)}")
    template = template.replace("orbit_number", f"{orbit:05d}")
    logger.info(f"scene: {cl.val(scene)}")
    template = template.replace("scene_ID", f"{scene:03d}")

    L2G_LSTE_filename = abspath(expanduser(L2G_LSTE_filename))

    logger.info(f"L2G_LSTE file: {cl.file(L2G_LSTE_filename)}")
    template = template.replace("L2G_LSTE_filename", L2G_LSTE_filename)

    if len(L2T_LSTE_filenames) == 0:
        raise ValueError(f"no L2T LSTE filenames given")

    logger.info(f"listing {len(L2T_LSTE_filenames)} L2T_LSTE files: ")

    L2T_LSTE_filenames_XML = "\n            ".join([
        f"<element>{abspath(expanduser(filename))}</element>"
        for filename
        in L2T_LSTE_filenames
    ])

    template = template.replace("<element>L2T_LSTE_filename1</element>", L2T_LSTE_filenames_XML)

    logger.info(f"working directory: {cl.dir(working_directory)}")
    template = template.replace("working_directory", working_directory)
    logger.info(f"L2T STARS sources directory: {cl.dir(L2T_STARS_sources_directory)}")
    template = template.replace("L2T_STARS_sources_directory", L2T_STARS_sources_directory)
    logger.info(f"L2T STARS indices directory: {cl.dir(L2T_STARS_indices_directory)}")
    template = template.replace("L2T_STARS_indices_directory", L2T_STARS_indices_directory)
    logger.info(f"L2T STARS model directory: {cl.dir(L2T_STARS_model_directory)}")
    template = template.replace("L2T_STARS_model_directory", L2T_STARS_model_directory)
    logger.info(f"executable: {cl.file(executable_filename)}")
    template = template.replace("executable_filename", executable_filename)
    logger.info(f"run-config: {cl.file(runconfig_filename)}")
    template = template.replace("runconfig_filename", runconfig_filename)
    logger.info(f"log: {cl.file(log_filename)}")
    template = template.replace("log_filename", log_filename)
    logger.info(f"build: {cl.val(build)}")
    template = template.replace("build_ID", build)
    logger.info(f"processing node: {cl.val(processing_node)}")
    template = template.replace("processing_node", processing_node)
    logger.info(f"production date/time: {cl.time(production_datetime)}")
    template = template.replace("production_datetime", production_datetime)
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

    return runconfig_filename


# FIXME still re-working from L3G_L4G to DOWNLOADER
class ECOv002DLConfig(ECOSTRESSRunConfig):
    def __init__(self, filename: str):
        try:
            logger.info(f"loading ECOv002_DL run-config: {cl.file(filename)}")
            runconfig = read_runconfig(filename)

            if "StaticAncillaryFileGroup" not in runconfig:
                raise MissingRunConfigValue(f"missing StaticAncillaryFileGroup in ECOv002_DL run-config: {filename}")

            if "ECOv002_DL_WORKING" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/ECOv002_DL_WORKING in ECOv002_DL run-config: {filename}")

            working_directory = abspath(runconfig["StaticAncillaryFileGroup"]["ECOv002_DL_WORKING"])
            logger.info(f"working directory: {cl.dir(working_directory)}")

            if "L2T_STARS_SOURCES" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2T_STARS_SOURCES in ECOv002_DL run-config: {filename}")

            L2T_STARS_sources_directory = abspath(runconfig["StaticAncillaryFileGroup"]["L2T_STARS_SOURCES"])
            logger.info(f"L2T STARS sources directory: {cl.dir(L2T_STARS_sources_directory)}")

            if "L2T_STARS_INDICES" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2T_STARS_INDICES in ECOv002_DL run-config: {filename}")

            L2T_STARS_indices_directory = abspath(runconfig["StaticAncillaryFileGroup"]["L2T_STARS_INDICES"])
            logger.info(f"L2T STARS indices directory: {cl.dir(L2T_STARS_indices_directory)}")

            if "L2T_STARS_MODEL" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2T_STARS_MODEL in ECOv002_DL run-config: {filename}")

            L2T_STARS_model_directory = abspath(runconfig["StaticAncillaryFileGroup"]["L2T_STARS_MODEL"])
            logger.info(f"L2T STARS model directory: {cl.dir(L2T_STARS_model_directory)}")

            if "ProductPathGroup" not in runconfig:
                raise MissingRunConfigValue(f"missing ProductPathGroup in ECOv002_DL run-config: {filename}")

            if "InputFileGroup" not in runconfig:
                raise MissingRunConfigValue(f"missing InputFileGroup in ECOv002_DL run-config: {filename}")

            if "L2G_LSTE" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(f"missing InputFileGroup/L2G_LSTE in ECOv002_DL run-config: {filename}")

            L2G_LSTE_filename = abspath(runconfig["InputFileGroup"]["L2G_LSTE"])
            logger.info(f"L2G_LSTE file: {cl.file(L2G_LSTE_filename)}")

            if "L2T_LSTE" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L2T_LSTE in ECOv002_DL run-config: {filename}")

            L2T_LSTE_filenames = runconfig["InputFileGroup"]["L2T_LSTE"]
            logger.info(f"reading {len(L2T_LSTE_filenames)} L2T_LSTE files")

            orbit = int(runconfig["Geometry"]["OrbitNumber"])
            logger.info(f"orbit: {cl.val(orbit)}")

            if "SceneId" not in runconfig["Geometry"]:
                raise MissingRunConfigValue(f"missing Geometry/SceneId in L2T_STARS run-config: {filename}")

            scene = int(runconfig["Geometry"]["SceneId"])
            logger.info(f"scene: {cl.val(scene)}")

            if "BuildID" not in runconfig["PrimaryExecutable"]:
                raise MissingRunConfigValue(
                    f"missing PrimaryExecutable/BuildID in L1_L2_RAD_LSTE run-config {filename}")

            build = str(runconfig["PrimaryExecutable"]["BuildID"])

            if "ProductCounter" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductCounter in L1_L2_RAD_LSTE run-config {filename}")

            product_counter = int(runconfig["ProductPathGroup"]["ProductCounter"])

            L2G_LSTE_granule = L2GLSTE(L2G_LSTE_filename)
            time_UTC = L2G_LSTE_granule.time_UTC

            timestamp = f"{time_UTC:%Y%m%dT%H%M%S}"

            PGE_name = "DOWNLOADER"
            PGE_version = ECOSTRESS.PGEVersion

            self.working_directory = working_directory
            self.L2G_LSTE_filename = L2G_LSTE_filename
            self.L2T_LSTE_filenames = L2T_LSTE_filenames
            self.L2T_STARS_sources_directory = L2T_STARS_sources_directory
            self.L2T_STARS_indices_directory = L2T_STARS_indices_directory
            self.L2T_STARS_model_directory = L2T_STARS_model_directory
            self.orbit = orbit
            self.scene = scene
            self.product_counter = product_counter
            self.time_UTC = time_UTC
            self.PGE_name = PGE_name
            self.PGE_version = PGE_version
        except MissingRunConfigValue as e:
            raise e
        except ECOSTRESSExitCodeException as e:
            raise e
        except Exception as e:
            logger.exception(e)
            raise UnableToParseRunConfig(f"unable to parse run-config file: {filename}")


def ECOv002_DL(runconfig_filename: str, tiles: List[str] = None) -> int:
    """
    ECOSTRESS Collection 2 Downloader PGE
    :param runconfig_filename: filename for XML run-config
    :return: exit code number
    """
    exit_code = SUCCESS_EXIT_CODE

    cl.configure()
    logger = logging.getLogger(__name__)

    try:
        logger.info(f"ECOSTRESS Collection 2 Downloader PGE ({cl.val(ECOSTRESS.PGEVersion)})")
        logger.info(f"run-config: {cl.file(runconfig_filename)}")
        runconfig = ECOv002DLConfig(runconfig_filename)
        working_directory = runconfig.working_directory
        logger.info(f"working directory: {cl.dir(working_directory)}")
        L2T_STARS_sources_directory = runconfig.L2T_STARS_sources_directory
        logger.info(f"L2T STARS sources directory: {cl.dir(L2T_STARS_sources_directory)}")
        L2T_STARS_indices_directory = runconfig.L2T_STARS_indices_directory
        logger.info(f"L2T STARS indices directory: {cl.dir(L2T_STARS_indices_directory)}")
        L2T_STARS_model_directory = runconfig.L2T_STARS_model_directory
        logger.info(f"L2T STARS model directory: {cl.dir(L2T_STARS_model_directory)}")
        L2G_LSTE_filename = runconfig.L2G_LSTE_filename
        logger.info(f"L2G LSTE file: {cl.file(L2G_LSTE_filename)}")

        orbit = runconfig.orbit
        logger.info(f"orbit: {cl.val(orbit)}")
        scene = runconfig.scene
        logger.info(f"scene: {cl.val(scene)}")

        logger.info(f"loading grid from L2G LSTE: {L2G_LSTE_filename}")
        L2G_LSTE_granule = L2GLSTE(L2G_LSTE_filename)
        geometry = L2G_LSTE_granule.grid
        time_UTC = L2G_LSTE_granule.time_UTC

        L2T_LSTE_filenames = runconfig.L2T_LSTE_filenames
        logger.info(
            f"processing {cl.val(len(L2T_LSTE_filenames))} tiles for orbit {cl.val(orbit)} scene {cl.val(scene)}")

        for L2T_LSTE_filename in L2T_LSTE_filenames:
            L2T_LSTE_granule = L2TLSTE(L2T_LSTE_filename)
            tile = L2T_LSTE_granule.tile
            tile_grid = SentinelTileGrid(target_resolution=70)

            if not tile_grid.land(tile):
                logger.warning(f"Sentinel tile {tile} is not on land")
                continue

            if tiles is not None and tile not in tiles:
                continue

            logger.info(f"L2T LSTE filename: {cl.file(L2T_LSTE_filename)}")
            logger.info(f"orbit: {cl.val(orbit)} scene: {cl.val(scene)} tile: {cl.val(tile)}")

            L2T_STARS_runconfig_filename = generate_L2T_STARS_runconfig(
                L2T_LSTE_filename=L2T_LSTE_filename,
                orbit=orbit,
                scene=scene,
                tile=tile,
                time_UTC=time_UTC,
                working_directory=working_directory,
                sources_directory=L2T_STARS_sources_directory,
                indices_directory=L2T_STARS_indices_directory,
                model_directory=L2T_STARS_model_directory
            )

            exit_code = L2T_STARS(
                runconfig_filename=L2T_STARS_runconfig_filename,
                sources_only=True
            )

            if exit_code == LAND_FILTER:
                logger.warning(f"Sentinel tile {tile} is not on land")
                continue

            if exit_code != 0:
                return exit_code

    except ECOSTRESSExitCodeException as exception:
        logger.exception(exception)
        exit_code = exception.exit_code

    return exit_code


def main(argv=sys.argv):
    if len(argv) == 1 or "--version" in argv:
        print(f"ECOSTRESS Collection 2 Downloader PGE ({ECOSTRESS.PGEVersion})")
        print(f"usage: ECOv002_DL RunConfig.xml")

        if "--version" in argv:
            return SUCCESS_EXIT_CODE
        else:
            return RUNCONFIG_FILENAME_NOT_SUPPLIED

    runconfig_filename = str(argv[1])
    exit_code = ECOv002_DL(runconfig_filename=runconfig_filename)
    logger.info(f"ECOSTRESS Collection 2 Downloader PGE exit code: {exit_code}")

    return exit_code


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
