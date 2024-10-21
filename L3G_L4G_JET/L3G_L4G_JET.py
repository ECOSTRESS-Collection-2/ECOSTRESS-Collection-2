import logging
import socket
import sys
from datetime import datetime
from os import makedirs
from os.path import join, abspath, dirname, expanduser, exists, basename, splitext
from shutil import which
from typing import List
from uuid import uuid4

import ECOSTRESS
import colored_logging as cl
from ECOSTRESS.L2_LSTE import L2GLSTE
from ECOSTRESS.L3_JET import L3GJET, L3GMET, L3GSEB, L3GSM
from ECOSTRESS.L4_ESI import L4GESI
from ECOSTRESS.L4_WUE import L4GWUE
from ECOSTRESS.exit_codes import SUCCESS_EXIT_CODE, ECOSTRESSExitCodeException, RUNCONFIG_FILENAME_NOT_SUPPLIED, \
    MissingRunConfigValue, UnableToParseRunConfig
from ECOSTRESS.runconfig import read_runconfig, ECOSTRESSRunConfig

with open(join(abspath(dirname(__file__)), "version.txt")) as f:
    version = f.read()

__version__ = version

L3G_L4G_JET_OUTPUT_DIRECTORY = "L3G_L4G_JET_output"
L3G_L4G_JET_TEMPLATE = join(abspath(dirname(__file__)), "L3G_L4G_JET.xml")
DEFAULT_BUILD = "0700"

L3G_SEB_SHORT_NAME = "ECO_L3G_SEB"
L3G_SEB_LONG_NAME = "ECOSTRESS Gridded Surface Energy Balance Instantaneous L3 Global 70 m"

L3G_SM_SHORT_NAME = "ECO_L3G_SM"
L3G_SM_LONG_NAME = "ECOSTRESS Gridded Downscaled Soil Moisture Instantaneous L3 Global 70 m"

L3G_MET_SHORT_NAME = "ECO_L3G_MET"
L3G_MET_LONG_NAME = "ECOSTRESS Gridded Downscaled Meteorology Instantaneous L3 Global 70 m"

L3G_JET_SHORT_NAME = "ECO_L3G_JET"
L3G_JET_LONG_NAME = "ECOSTRESS Gridded Evapotranspiration Instantaneous and Daytime L3 Global 70 m"

L4G_ESI_SHORT_NAME = "ECO_L4G_ESI"
L4G_ESI_LONG_NAME = "ECOSTRESS Gridded Evaporative Stress Index Instantaneous L4 Global 70 m"

L4G_WUE_SHORT_NAME = "ECO_L4G_WUE"
L4G_WUE_LONG_NAME = "ECOSTRESS Gridded Water Use Efficiency Instantaneous L4 Global 70 m"

logger = logging.getLogger(__name__)


def generate_L3G_L4G_JET_runconfig(
        L2G_LSTE_filename: str,
        L3T_JET_filenames: List[str],
        L3T_MET_filenames: List[str],
        L3T_SEB_filenames: List[str],
        L3T_SM_filenames: List[str],
        L4T_ESI_filenames: List[str],
        L4T_WUE_filenames: List[str],
        orbit: int = None,
        scene: int = None,
        working_directory: str = None,
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
        template_filename = L3G_L4G_JET_TEMPLATE

    template_filename = abspath(expanduser(template_filename))

    run_ID = f"ECOv002_L3G_L4G_JET_{orbit:05d}_{scene:05d}"

    if runconfig_filename is None:
        runconfig_filename = join(working_directory, "runconfig", f"{run_ID}.xml")

    runconfig_filename = abspath(expanduser(runconfig_filename))

    if working_directory is None:
        working_directory = run_ID

    working_directory = abspath(expanduser(working_directory))

    if executable_filename is None:
        executable_filename = which("L3G_L4G_JET")

    if executable_filename is None:
        executable_filename = "L3G_L4G_JET"

    if output_directory is None:
        output_directory = join(working_directory, L3G_L4G_JET_OUTPUT_DIRECTORY)

    output_directory = abspath(expanduser(output_directory))

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
    logger.info(f"loading L3G_L4G_JET template: {cl.file(template_filename)}")

    with open(template_filename, "r") as file:
        template = file.read()

    logger.info(f"orbit: {cl.val(orbit)}")
    template = template.replace("orbit_number", f"{orbit:05d}")
    logger.info(f"scene: {cl.val(scene)}")
    template = template.replace("scene_ID", f"{scene:03d}")

    L2G_LSTE_filename = abspath(expanduser(L2G_LSTE_filename))

    logger.info(f"L2G_LSTE file: {cl.file(L2G_LSTE_filename)}")
    template = template.replace("L2G_LSTE_filename", L2G_LSTE_filename)

    if len(L3T_JET_filenames) == 0:
        raise ValueError(f"no L3T JET filenames given")

    logger.info(f"listing {len(L3T_JET_filenames)} L3T_JET files: ")

    L3T_JET_filenames_XML = "\n            ".join([
        f"<element>{abspath(expanduser(filename))}</element>"
        for filename
        in L3T_JET_filenames
    ])

    template = template.replace("<element>L3T_JET_filename1</element>", L3T_JET_filenames_XML)

    if len(L3T_MET_filenames) == 0:
        raise ValueError("no L3T MET filenames given")

    logger.info(f"listing {len(L3T_MET_filenames)} L3T_MET files: ")

    L3T_MET_filenames_XML = "\n            ".join([
        f"<element>{abspath(expanduser(filename))}</element>"
        for filename
        in L3T_MET_filenames
    ])

    template = template.replace("<element>L3T_MET_filename1</element>", L3T_MET_filenames_XML)

    if len(L3T_SEB_filenames) == 0:
        raise ValueError(f"no L3T SEB filenames given")

    logger.info(f"listing {len(L3T_SEB_filenames)} L3T_SEB files: ")

    L3T_SEB_filenames_XML = "\n            ".join([
        f"<element>{abspath(expanduser(filename))}</element>"
        for filename
        in L3T_SEB_filenames
    ])

    if len(L3T_SM_filenames) == 0:
        raise ValueError(f"no L3T SM filenames given")

    logger.info(f"listing {len(L3T_SM_filenames)} L3T_SM files: ")

    L3T_SM_filenames_XML = "\n            ".join([
        f"<element>{abspath(expanduser(filename))}</element>"
        for filename
        in L3T_SM_filenames
    ])

    template = template.replace("<element>L3T_SM_filename1</element>", L3T_SM_filenames_XML)

    template = template.replace("<element>L3T_SEB_filename1</element>", L3T_SEB_filenames_XML)

    if len(L4T_ESI_filenames) == 0:
        raise ValueError(f"no L4T ESI filenames given")

    logger.info(f"listing {len(L4T_ESI_filenames)} L4T_ESI files: ")

    L4T_ESI_filenames_XML = "\n            ".join([
        f"<element>{abspath(expanduser(filename))}</element>"
        for filename
        in L4T_ESI_filenames
    ])

    template = template.replace("<element>L4T_ESI_filename1</element>", L4T_ESI_filenames_XML)

    if len(L4T_WUE_filenames) == 0:
        raise ValueError(f"no L4T WUE filenames given")

    logger.info(f"listing {len(L4T_WUE_filenames)} L4T_WUE files: ")

    L4T_WUE_filenames_XML = "\n            ".join([
        f"<element>{abspath(expanduser(filename))}</element>"
        for filename
        in L4T_WUE_filenames
    ])

    template = template.replace("<element>L4T_WUE_filename1</element>", L4T_WUE_filenames_XML)

    logger.info(f"working directory: {cl.dir(working_directory)}")
    template = template.replace("working_directory", working_directory)
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

class L3GL4GJETConfig(ECOSTRESSRunConfig):
    def __init__(self, filename: str):
        try:
            logger.info(f"loading L3G_L4G_JET run-config: {cl.file(filename)}")
            runconfig = read_runconfig(filename)

            if "StaticAncillaryFileGroup" not in runconfig:
                raise MissingRunConfigValue(f"missing StaticAncillaryFileGroup in L3G_L4G_JET run-config: {filename}")

            if "L3G_L4G_JET_WORKING" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L3G_L4G_JET_WORKING in L3G_L4G_JET run-config: {filename}")

            working_directory = abspath(runconfig["StaticAncillaryFileGroup"]["L3G_L4G_JET_WORKING"])
            logger.info(f"working directory: {cl.dir(working_directory)}")

            if "ProductPathGroup" not in runconfig:
                raise MissingRunConfigValue(f"missing ProductPathGroup in L3G_L4G_JET run-config: {filename}")

            if "ProductPath" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductPath in L3G_L4G_JET run-config: {filename}")

            output_directory = abspath(runconfig["ProductPathGroup"]["ProductPath"])
            logger.info(f"output directory: {cl.dir(output_directory)}")

            if "InputFileGroup" not in runconfig:
                raise MissingRunConfigValue(f"missing InputFileGroup in L3G_L4G_JET run-config: {filename}")

            if "L2G_LSTE" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(f"missing InputFileGroup/L2G_LSTE in L3G_L4G_JET run-config: {filename}")

            L2G_LSTE_filename = abspath(runconfig["InputFileGroup"]["L2G_LSTE"])
            logger.info(f"L2G_LSTE file:{cl.file(L2G_LSTE_filename)}")

            if "L3T_JET" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L3T_JET in L3G_L4G_JET run-config: {filename}")

            L3T_JET_filenames = runconfig["InputFileGroup"]["L3T_JET"]
            logger.info(f"reading {len(L3T_JET_filenames)} L3T_JET files")

            if "L3T_MET" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L3T_MET in L3G_L4G_JET run-config: {filename}")

            L3T_MET_filenames = runconfig["InputFileGroup"]["L3T_MET"]
            logger.info(f"reading {len(L3T_MET_filenames)} L3T_MET files")

            if "L3T_SEB" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L3T_SEB in L3G_L4G run-config: {filename}")

            L3T_SEB_filenames = runconfig["InputFileGroup"]["L3T_SEB"]
            logger.info(f"reading {len(L3T_SEB_filenames)} L3T_SEB files")

            if "L3T_SM" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L3T_SM in L3G_L4G_JET run-config: {filename}")

            L3T_SM_filenames = runconfig["InputFileGroup"]["L3T_SM"]
            logger.info(f"reading {len(L3T_SM_filenames)} L3T_SM files")

            if "L4T_ESI" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L4T_ESI in L3G_L4G_JET run-config: {filename}")

            L4T_ESI_filenames = runconfig["InputFileGroup"]["L4T_ESI"]
            logger.info(f"reading {len(L4T_ESI_filenames)} L4T_ESI files")

            if "L4T_WUE" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(f"missing InputFileGroup/L4T_WUE in L3G_L4G_JET run-config: {filename}")

            L4T_WUE_filenames = runconfig["InputFileGroup"]["L4T_WUE"]
            logger.info(f"reading {len(L4T_WUE_filenames)} L4T_WUE files")

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

            L3G_JET_granule_ID = f"ECOv002_L3G_JET_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L3G_JET_filename = join(output_directory, f"{L3G_JET_granule_ID}.h5")

            L3G_MET_granule_ID = f"ECOv002_L3G_MET_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L3G_MET_filename = join(output_directory, f"{L3G_MET_granule_ID}.h5")

            L3G_SEB_granule_ID = f"ECOv002_L3G_SEB_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L3G_SEB_filename = join(output_directory, f"{L3G_SEB_granule_ID}.h5")

            L3G_SM_granule_ID = f"ECOv002_L3G_SM_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L3G_SM_filename = join(output_directory, f"{L3G_SM_granule_ID}.h5")

            L4G_ESI_granule_ID = f"ECOv002_L4G_ESI_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L4G_ESI_filename = join(output_directory, f"{L4G_ESI_granule_ID}.h5")

            L4G_WUE_granule_ID = f"ECOv002_L4G_WUE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L4G_WUE_filename = join(output_directory, f"{L4G_WUE_granule_ID}.h5")

            PGE_name = "L3G_L4G_JET"
            PGE_version = ECOSTRESS.PGEVersion

            self.working_directory = working_directory
            self.output_directory = output_directory
            self.L2G_LSTE_filename = L2G_LSTE_filename

            self.L3T_JET_filenames = L3T_JET_filenames
            self.L3G_JET_granule_ID = L3G_JET_granule_ID
            self.L3G_JET_filename = L3G_JET_filename

            self.L3T_MET_filenames = L3T_MET_filenames
            self.L3G_MET_granule_ID = L3G_MET_granule_ID
            self.L3G_MET_filename = L3G_MET_filename

            self.L3T_SEB_filenames = L3T_SEB_filenames
            self.L3G_SEB_granule_ID = L3G_SEB_granule_ID
            self.L3G_SEB_filename = L3G_SEB_filename

            self.L3T_SM_filenames = L3T_SM_filenames
            self.L3G_SM_granule_ID = L3G_SM_granule_ID
            self.L3G_SM_filename = L3G_SM_filename

            self.L4T_ESI_filenames = L4T_ESI_filenames
            self.L4G_ESI_granule_ID = L4G_ESI_granule_ID
            self.L4G_ESI_filename = L4G_ESI_filename

            self.L4T_WUE_filenames = L4T_WUE_filenames
            self.L4G_WUE_granule_ID = L4G_WUE_granule_ID
            self.L4G_WUE_filename = L4G_WUE_filename

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


def L3G_L4G_JET(runconfig_filename: str, geotiff_diagnostics: bool = False) -> int:
    """
    ECOSTRESS Collection 2 L3G L4G JET PGE
    :param runconfig_filename: filename for XML run-config
    :param log_filename: filename for logger output
    :return: exit code number
    """
    exit_code = SUCCESS_EXIT_CODE

    cl.configure()
    logger = logging.getLogger(__name__)

    try:
        logger.info(f"L3G L4G JET PGE ({cl.val(ECOSTRESS.PGEVersion)})")
        logger.info(f"run-config: {cl.file(runconfig_filename)}")
        runconfig = L3GL4GJETConfig(runconfig_filename)
        working_directory = runconfig.working_directory
        L2G_LSTE_filename = runconfig.L2G_LSTE_filename
        L3T_JET_filenames = runconfig.L3T_JET_filenames
        L3G_JET_filename = runconfig.L3G_JET_filename
        L3T_MET_filenames = runconfig.L3T_MET_filenames
        L3G_MET_filename = runconfig.L3G_MET_filename
        L3T_SEB_filenames = runconfig.L3T_SEB_filenames
        L3G_SEB_filename = runconfig.L3G_SEB_filename
        L3T_SM_filenames = runconfig.L3T_SM_filenames
        L3G_SM_filename = runconfig.L3G_SM_filename
        L4T_ESI_filenames = runconfig.L4T_ESI_filenames
        L4G_ESI_filename = runconfig.L4G_ESI_filename
        L4T_WUE_filenames = runconfig.L4T_WUE_filenames
        L4G_WUE_filename = runconfig.L4G_WUE_filename
        orbit = runconfig.orbit
        scene = runconfig.scene

        logger.info(f"loading grid from L2G LSTE: {L2G_LSTE_filename}")
        L2G_LSTE_granule = L2GLSTE(L2G_LSTE_filename)
        geometry = L2G_LSTE_granule.grid

        standard_metadata = {
            "SISName": "Level 3/4 JET Product Specification Document",
            "SISVersion": "Preliminary",
            "PGEName": "L3G_L4G_JET"
        }

        L3G_JET_browse_filename = L3G_JET_filename.replace(".h5", ".png")

        if exists(L3G_JET_filename) and exists(L3G_JET_browse_filename):
            L3G_JET_granule = L3GJET(L3G_JET_filename)
        else:
            standard_metadata["ProcessingLevelID"] = "L3G"
            standard_metadata["ProcessingLevelDescription"] = "Level 3 Gridded Evapotranspiration Ensemble"

            short_name = L3G_JET_SHORT_NAME
            logger.info(f"L3G JET short name: {cl.name(short_name)}")
            standard_metadata["ShortName"] = short_name

            long_name = L3G_JET_LONG_NAME
            logger.info(f"L3G JET long name: {cl.name(long_name)}")
            standard_metadata["LongName"] = long_name

            logger.info(f"generating L3G JET granule: {cl.file(L3G_JET_filename)}")
            L3G_JET_granule = L3GJET.from_tiles(
                filename=L3G_JET_filename,
                tile_filenames=L3T_JET_filenames,
                gridded_source_granule=L2G_LSTE_granule,
                standard_metadata_additional=standard_metadata,
                geotiff_diagnostics=geotiff_diagnostics
            )

            logger.info(f"generating L3G JET browse image: {cl.file(L3G_JET_browse_filename)}")
            L3G_JET_granule.write_browse_image(PNG_filename=L3G_JET_browse_filename)

        L3G_MET_browse_filename = L3G_MET_filename.replace(".h5", ".png")

        if exists(L3G_MET_filename) and exists(L3G_MET_browse_filename):
            L3G_MET_granule = L3GMET(L3G_MET_filename)
        else:
            standard_metadata["ProcessingLevelID"] = "L3G"
            standard_metadata["ProcessingLevelDescription"] = "Level 3 Gridded Meteorology"

            short_name = L3G_MET_SHORT_NAME
            logger.info(f"L3G MET short name: {cl.name(short_name)}")
            standard_metadata["ShortName"] = short_name

            long_name = L3G_MET_LONG_NAME
            logger.info(f"L3G MET long name: {cl.name(long_name)}")
            standard_metadata["LongName"] = long_name

            logger.info(f"generating L3G MET granule: {cl.file(L3G_MET_filename)}")
            L3G_MET_granule = L3GMET.from_tiles(
                filename=L3G_MET_filename,
                tile_filenames=L3T_MET_filenames,
                gridded_source_granule=L2G_LSTE_granule,
                standard_metadata_additional=standard_metadata,
                geotiff_diagnostics=geotiff_diagnostics
            )

            logger.info(f"generating L3G MET browse image: {cl.file(L3G_MET_browse_filename)}")
            L3G_MET_granule.write_browse_image(PNG_filename=L3G_MET_browse_filename)

        L3G_SEB_browse_filename = L3G_SEB_filename.replace(".h5", ".png")

        if exists(L3G_SEB_filename) and exists(L3G_SEB_browse_filename):
            L3G_SEB_granule = L3GSEB(L3G_SEB_filename)
        else:
            standard_metadata["ProcessingLevelID"] = "L3G"
            standard_metadata["ProcessingLevelDescription"] = "Level 3 Gridded Surface Energy Balance"

            short_name = L3G_SEB_SHORT_NAME
            logger.info(f"L3G SEB short name: {cl.name(short_name)}")
            standard_metadata["ShortName"] = short_name

            long_name = L3G_SEB_LONG_NAME
            logger.info(f"L3G SEB long name: {cl.name(long_name)}")
            standard_metadata["LongName"] = long_name

            logger.info(f"generating L3G SEB granule: {cl.file(L3G_SEB_filename)}")
            L3G_SEB_granule = L3GSEB.from_tiles(
                filename=L3G_SEB_filename,
                tile_filenames=L3T_SEB_filenames,
                gridded_source_granule=L2G_LSTE_granule,
                standard_metadata_additional=standard_metadata,
                geotiff_diagnostics=geotiff_diagnostics
            )

            logger.info(f"generating L3G SEB browse image: {cl.file(L3G_SEB_browse_filename)}")
            L3G_SEB_granule.write_browse_image(PNG_filename=L3G_SEB_browse_filename)

        L3G_SM_browse_filename = L3G_SM_filename.replace(".h5", ".png")

        if exists(L3G_SM_filename) and exists(L3G_SM_browse_filename):
            L3G_SM_granule = L3GSM(L3G_SM_filename)
        else:
            standard_metadata["ProcessingLevelID"] = "L3G"
            standard_metadata["ProcessingLevelDescription"] = "Level 3 Gridded Soil Moisture"

            short_name = L3G_SM_SHORT_NAME
            logger.info(f"L3G SM short name: {cl.name(short_name)}")
            standard_metadata["ShortName"] = short_name

            long_name = L3G_SM_LONG_NAME
            logger.info(f"L3G SM long name: {cl.name(long_name)}")
            standard_metadata["LongName"] = long_name

            logger.info(f"generating L3G SM granule: {cl.file(L3G_SM_filename)}")
            L3G_SM_granule = L3GSM.from_tiles(
                filename=L3G_SM_filename,
                tile_filenames=L3T_SM_filenames,
                gridded_source_granule=L2G_LSTE_granule,
                standard_metadata_additional=standard_metadata,
                geotiff_diagnostics=geotiff_diagnostics
            )

            logger.info(f"generating L3G SM browse image: {cl.file(L3G_SM_browse_filename)}")
            L3G_SM_granule.write_browse_image(PNG_filename=L3G_SM_browse_filename)

        L4G_ESI_browse_filename = L4G_ESI_filename.replace(".h5", ".png")

        if exists(L4G_ESI_filename) and exists(L4G_ESI_browse_filename):
            L4G_ESI_granule = L4GESI(L4G_ESI_filename)
        else:
            standard_metadata["ProcessingLevelID"] = "L4G"
            standard_metadata["ProcessingLevelDescription"] = "Level 4 Gridded Evaporative Stress Index"

            short_name = L4G_ESI_SHORT_NAME
            logger.info(f"L4G ESI short name: {cl.name(short_name)}")
            standard_metadata["ShortName"] = short_name

            long_name = L4G_ESI_LONG_NAME
            logger.info(f"L4G ESI long name: {cl.name(long_name)}")
            standard_metadata["LongName"] = long_name

            logger.info(f"generating L4G ESI granule: {cl.file(L4G_ESI_filename)}")
            L4G_ESI_granule = L4GESI.from_tiles(
                filename=L4G_ESI_filename,
                tile_filenames=L4T_ESI_filenames,
                gridded_source_granule=L2G_LSTE_granule,
                standard_metadata_additional=standard_metadata,
                geotiff_diagnostics=geotiff_diagnostics
            )

            logger.info(f"generating L4G ESI browse image: {cl.file(L4G_ESI_browse_filename)}")
            L4G_ESI_granule.write_browse_image(PNG_filename=L4G_ESI_browse_filename)

        L4G_WUE_browse_filename = L4G_WUE_filename.replace(".h5", ".png")

        if exists(L4G_WUE_filename) and exists(L4G_ESI_browse_filename):
            L4G_WUE_granule = L4GWUE(L4G_WUE_filename)
        else:
            standard_metadata["ProcessingLevelID"] = "L4G"
            standard_metadata["ProcessingLevelDescription"] = "Level 4 Gridded Water Use Efficiency"

            short_name = L4G_WUE_SHORT_NAME
            logger.info(f"L4G WUE short name: {cl.name(short_name)}")
            standard_metadata["ShortName"] = short_name

            long_name = L4G_WUE_LONG_NAME
            logger.info(f"L4G WUE long name: {cl.name(long_name)}")
            standard_metadata["LongName"] = long_name

            logger.info(f"generating L4G WUE granule: {cl.file(L4G_WUE_filename)}")
            L4G_WUE_granule = L4GWUE.from_tiles(
                filename=L4G_WUE_filename,
                tile_filenames=L4T_WUE_filenames,
                gridded_source_granule=L2G_LSTE_granule,
                standard_metadata_additional=standard_metadata,
                geotiff_diagnostics=geotiff_diagnostics
            )

            logger.info(f"generating L4G WUE browse image: {cl.file(L4G_WUE_browse_filename)}")
            L4G_WUE_granule.write_browse_image(PNG_filename=L4G_WUE_browse_filename)

    except ECOSTRESSExitCodeException as exception:
        logger.exception(exception)
        exit_code = exception.exit_code

    return exit_code


def main(argv=sys.argv):
    if len(argv) == 1 or "--version" in argv:
        print(f"L3G_L4G_JET PGE ({ECOSTRESS.PGEVersion})")
        print(f"usage: L3G_L4G_JET RunConfig.xml")

        if "--version" in argv:
            return SUCCESS_EXIT_CODE
        else:
            return RUNCONFIG_FILENAME_NOT_SUPPLIED

    runconfig_filename = str(argv[1])
    exit_code = L3G_L4G_JET(runconfig_filename=runconfig_filename)
    logger.info(f"L3G_L4G_JET exit code: {exit_code}")

    return exit_code


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
