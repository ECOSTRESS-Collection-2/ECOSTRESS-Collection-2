import logging
import socket
import sys
from datetime import datetime
from glob import glob
from os import makedirs
from os.path import join, abspath, dirname, expanduser, exists, basename, splitext
from shutil import which
from typing import List
from uuid import uuid4

from dateutil import parser

import ECOSTRESS
import colored_logging as cl
from ECOSTRESS.L1_RAD import L1CGRAD, L1BRAD
from ECOSTRESS.L2G_CLOUD import L2GCLOUD
from ECOSTRESS.L2_CLOUD import L2CLOUD
from ECOSTRESS.L2_LSTE import L2LSTE, L2GLSTE
from ECOSTRESS.exit_codes import SUCCESS_EXIT_CODE, ECOSTRESSExitCodeException, RUNCONFIG_FILENAME_NOT_SUPPLIED, \
    MissingRunConfigValue, UnableToParseRunConfig, LandFilter
from ECOSTRESS.find_ECOSTRESS_C1_scene import find_ECOSTRESS_C1_scene
from ECOSTRESS.runconfig import ECOSTRESSRunConfig
from ECOSTRESS.scan_resampling import generate_scan_kd_trees, clip_tails
from rasters import KDTree, RasterGrid
from timer import Timer

with open(join(abspath(dirname(__file__)), "version.txt")) as f:
    version = f.read()

__version__ = version

PGE_NAME = "L1_L2_RAD_LSTE"
L1_L2_RAD_LSTE_TEMPLATE = join(abspath(dirname(__file__)), "L1_L2_RAD_LSTE.xml")
DEFAULT_BUILD = "0700"
DEFAULT_WORKING_DIRECTORY = "."
DEFAULT_OUTPUT_DIRECTORY = "L1_L2_RAD_LSTE_output"
STRIP_CONSOLE = False
CELL_SIZE_DEGREES = 0.0007
CELL_SIZE_METERS = 70
SEARCH_RADIUS_METERS = 100
PROJECTION_SYSTEM = "global_geographic"
OVERLAP_STRATEGY = "checkerboard"

logger = logging.getLogger(__name__)


def generate_L1_L2_RAD_LSTE_runconfig(
        L2_LSTE_filename: str,
        L2_CLOUD_filename: str = None,
        L1B_GEO_filename: str = None,
        L1B_RAD_filename: str = None,
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
    L2_LSTE_filename = abspath(expanduser(L2_LSTE_filename))

    if not exists(L2_LSTE_filename):
        raise IOError(f"L2 LSTE file not found: {L2_LSTE_filename}")

    logger.info(f"L2 LSTE file: {cl.file(L2_LSTE_filename)}")
    source_granule_ID = splitext(basename(L2_LSTE_filename))[0]
    logger.info(f"source granule ID: {cl.name(source_granule_ID)}")

    if orbit is None:
        orbit = int(source_granule_ID.split("_")[-5])

    logger.info(f"orbit: {cl.val(orbit)}")

    if scene is None:
        scene = int(source_granule_ID.split("_")[-4])

    logger.info(f"scene: {cl.val(scene)}")

    if build is None:
        build = source_granule_ID.split("_")[-2]

    if L2_CLOUD_filename is None:
        directory = abspath(expanduser(dirname(L2_LSTE_filename)))
        pattern = join(directory, f"*_L2_CLOUD_{orbit:05d}_{scene:03d}_*_{build}_*.h5")
        logger.info(f"searching for L2 CLOUD: {cl.val(pattern)}")
        candidates = sorted(glob(pattern))

        if len(candidates) == 0:
            pattern = join(directory, f"*_L2_CLOUD_{orbit:05d}_{scene:03d}_*.h5")
            logger.info(f"searching for L2 CLOUD: {cl.val(pattern)}")
            candidates = sorted(glob(pattern))

        if len(candidates) == 0:
            raise ValueError("no L2 CLOUD filename given or found")

        L2_CLOUD_filename = candidates[-1]

    logger.info(f"L2 CLOUD file: {cl.file(L2_CLOUD_filename)}")

    if L1B_GEO_filename is None:
        directory = abspath(expanduser(dirname(L2_LSTE_filename)))
        pattern = join(directory, f"*_L1B_GEO_{orbit:05d}_{scene:03d}_*_{build}_*.h5")
        logger.info(f"searching for L1B GEO: {cl.val(pattern)}")
        candidates = sorted(glob(pattern))

        if len(candidates) == 0:
            pattern = join(directory, f"*_L1B_GEO_{orbit:05d}_{scene:03d}_*.h5")
            logger.info(f"searching for L1B GEO: {cl.val(pattern)}")
            candidates = sorted(glob(pattern))

        if len(candidates) == 0:
            raise ValueError("no L1B GEO filename given or found")

        L1B_GEO_filename = candidates[-1]

    logger.info(f"L1B GEO file: {cl.file(L1B_GEO_filename)}")

    if L1B_RAD_filename is None:
        directory = abspath(expanduser(dirname(L2_LSTE_filename)))
        pattern = join(directory, f"*_L1B_RAD_{orbit:05d}_{scene:03d}_*_{build}_*.h5")
        logger.info(f"searching for L1B RAD: {cl.val(pattern)}")
        candidates = sorted(glob(pattern))

        if len(candidates) == 0:
            pattern = join(directory, f"*_L1B_RAD_{orbit:05d}_{scene:03d}_*.h5")
            logger.info(f"searching for L1B RAD: {cl.val(pattern)}")
            candidates = sorted(glob(pattern))

        if len(candidates) == 0:
            raise ValueError("no L1B RAD filename given or found")

        L1B_RAD_filename = candidates[-1]

    logger.info(f"L1B RAD file: {cl.file(L1B_RAD_filename)}")

    if template_filename is None:
        template_filename = L1_L2_RAD_LSTE_TEMPLATE

    template_filename = abspath(expanduser(template_filename))

    if executable_filename is None:
        executable_filename = which("L1_L2_RAD_LSTE")

    if executable_filename is None:
        executable_filename = "L1_L2_RAD_LSTE"

    if output_directory is None:
        output_directory = join(working_directory, DEFAULT_OUTPUT_DIRECTORY)

    output_directory = abspath(expanduser(output_directory))

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

    L2_LSTE_filename = abspath(expanduser(L2_LSTE_filename))
    L2_CLOUD_filename = abspath(expanduser(L2_CLOUD_filename))
    L1B_GEO_filename = abspath(expanduser(L1B_GEO_filename))

    L2_LSTE_granule = L2LSTE(
        L2_LSTE_filename=L2_LSTE_filename,
        L2_CLOUD_filename=L2_CLOUD_filename,
        L1B_GEO_filename=L1B_GEO_filename
    )

    time_UTC = L2_LSTE_granule.time_UTC

    timestamp = f"{time_UTC:%Y%m%dT%H%M%S}"
    granule_ID = f"ECOv002_L1_L2_RAD_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"

    if runconfig_filename is None:
        runconfig_filename = join(working_directory, "runconfig", f"{granule_ID}.xml")

    runconfig_filename = abspath(expanduser(runconfig_filename))

    if exists(runconfig_filename):
        return runconfig_filename

    if log_filename is None:
        log_filename = join(working_directory, "log", f"{granule_ID}.log")

    log_filename = abspath(expanduser(log_filename))

    if working_directory is None:
        working_directory = granule_ID

    working_directory = abspath(expanduser(working_directory))

    logger.info(f"generating run-config for orbit {cl.val(orbit)} scene {cl.val(scene)}")
    logger.info(f"loading L1_L2_RAD_LSTE template: {cl.file(template_filename)}")

    with open(template_filename, "r") as file:
        template = file.read()

    logger.info(f"orbit: {cl.val(orbit)}")
    template = template.replace("orbit_number", f"{orbit:05d}")
    logger.info(f"scene: {cl.val(scene)}")
    template = template.replace("scene_ID", f"{scene:03d}")
    logger.info(f"L2_LSTE file: {cl.file(L2_LSTE_filename)}")
    template = template.replace("L2_LSTE_filename", L2_LSTE_filename)
    logger.info(f"L2_CLOUD file: {cl.file(L2_CLOUD_filename)}")
    template = template.replace("L2_CLOUD_filename", L2_CLOUD_filename)
    logger.info(f"L1B_GEO file: {cl.file(L1B_GEO_filename)}")
    template = template.replace("L1B_GEO_filename", L1B_GEO_filename)
    logger.info(f"L1B_RAD file: {cl.file(L1B_RAD_filename)}")
    template = template.replace("L1B_RAD_filename", L1B_RAD_filename)
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


def L1_L2_RAD_LSTE_runconfig_from_C1(orbit: int, scene: int, runconfig_filename: str) -> str:
    filenames = find_ECOSTRESS_C1_scene(
        orbit=orbit,
        scene=scene
    )

    L2_LSTE_filename = filenames["L2_LSTE"]
    L2_CLOUD_filename = filenames["L2_CLOUD"]
    L1B_GEO_filename = filenames["L1B_GEO"]

    runconfig = generate_L1_L2_RAD_LSTE_runconfig(
        orbit=orbit,
        scene=scene,
        L2_LSTE_filename=L2_LSTE_filename,
        L2_CLOUD_filename=L2_CLOUD_filename,
        L1B_GEO_filename=L1B_GEO_filename,
        runconfig_filename=runconfig_filename
    )

    return runconfig


class L2GL2TRADLSTEConfig(ECOSTRESSRunConfig):
    def __init__(self, filename: str):
        logger.info(f"loading L1_L2_RAD_LSTE run-config: {cl.file(filename)}")
        runconfig = self.read_runconfig(filename)

        # print(JSON_highlight(runconfig))

        try:
            if "StaticAncillaryFileGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup in L1_L2_RAD_LSTE run-config: {filename}")

            if "L2G_L2T_WORKING" not in runconfig["StaticAncillaryFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup/L2G_L2T_WORKING in L1_L2_RAD_LSTE run-config: {filename}")  # TODO exit code

            working_directory = abspath(runconfig["StaticAncillaryFileGroup"]["L2G_L2T_WORKING"])
            logger.info(f"working directory: {cl.dir(working_directory)}")

            if "ProductPathGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup in L1_L2_RAD_LSTE run-config: {filename}")

            if "ProductPath" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductPath in L1_L2_RAD_LSTE run-config: {filename}")

            output_directory = abspath(runconfig["ProductPathGroup"]["ProductPath"])
            logger.info(f"output directory: {cl.dir(output_directory)}")

            if "InputFileGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup in L1_L2_RAD_LSTE run-config: {filename}")

            if "L2_LSTE" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L2_LSTE in L1_L2_RAD_LSTE run-config: {filename}")

            L2_LSTE_filename = abspath(runconfig["InputFileGroup"]["L2_LSTE"])
            logger.info(f"L2_LSTE file: {cl.file(L2_LSTE_filename)}")

            if "L2_CLOUD" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L2_CLOUD in L1_L2_RAD_LSTE run-config: {filename}")

            L2_CLOUD_filename = abspath(runconfig["InputFileGroup"]["L2_CLOUD"])
            logger.info(f"L2_CLOUD file: {cl.file(L2_CLOUD_filename)}")

            if "L1B_GEO" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L1B_GEO in L1_L2_RAD_LSTE run-config: {filename}")

            L1B_GEO_filename = abspath(runconfig["InputFileGroup"]["L1B_GEO"])
            logger.info(f"L1B_GEO file: {cl.file(L1B_GEO_filename)}")

            if "L1B_RAD" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L1B_RAD in L1_L2_RAD_LSTE run-config: {filename}")

            L1B_RAD_filename = abspath(runconfig["InputFileGroup"]["L1B_RAD"])
            logger.info(f"L1B_RAD file: {cl.file(L1B_RAD_filename)}")

            orbit = int(runconfig["Geometry"]["OrbitNumber"])
            logger.info(f"orbit: {cl.val(orbit)}")

            if "SceneId" not in runconfig["Geometry"]:
                raise MissingRunConfigValue(
                    f"missing Geometry/SceneId in L1_L2_RAD_LSTE run-config: {filename}")

            scene = int(runconfig["Geometry"]["SceneId"])
            logger.info(f"scene: {cl.val(scene)}")

            if "ProductionDateTime" not in runconfig["JobIdentification"]:
                raise MissingRunConfigValue(
                    f"missing JobIdentification/ProductionDateTime in L1_L2_RAD_LSTE run-config {filename}")

            production_datetime = parser.parse(runconfig["JobIdentification"]["ProductionDateTime"])
            logger.info(f"production time: {cl.time(production_datetime)}")

            if "BuildID" not in runconfig["PrimaryExecutable"]:
                raise MissingRunConfigValue(
                    f"missing PrimaryExecutable/BuildID in L1_L2_RAD_LSTE run-config {filename}")

            build = str(runconfig["PrimaryExecutable"]["BuildID"])

            if "ProductCounter" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductCounter in L1_L2_RAD_LSTE run-config {filename}")

            product_counter = int(runconfig["ProductPathGroup"]["ProductCounter"])

            L2_LSTE_granule = L2LSTE(
                L2_LSTE_filename=L2_LSTE_filename,
                L2_CLOUD_filename=L2_CLOUD_filename,
                L1B_GEO_filename=L1B_GEO_filename
            )

            time_UTC = L2_LSTE_granule.time_UTC
            land_percent = L2_LSTE_granule.land_percent

            if land_percent == 0:
                raise LandFilter(f"skipping ocean scene with OverAllLandFraction value of {land_percent}")

            timestamp = f"{time_UTC:%Y%m%dT%H%M%S}"
            granule_ID = f"ECOv002_L1_L2_RAD_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L1CG_RAD_granule_ID = f"ECOv002_L1CG_RAD_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L2G_LSTE_granule_ID = f"ECOv002_L2G_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L2G_CLOUD_granule_ID = f"ECOv002_L2G_CLOUD_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            log_filename = abspath(expanduser(join(working_directory, "log", f"{granule_ID}.log")))
            L1CG_RAD_filename = join(output_directory, f"{L1CG_RAD_granule_ID}.h5")
            L2G_LSTE_filename = join(output_directory, f"{L2G_LSTE_granule_ID}.h5")
            L2G_CLOUD_filename = join(output_directory, f"{L2G_CLOUD_granule_ID}.h5")
            PGE_name = PGE_NAME
            PGE_version = ECOSTRESS.PGEVersion

            self.working_directory = working_directory
            self.log_filename = log_filename
            self.output_directory = output_directory
            self.L2_LSTE_filename = L2_LSTE_filename
            self.L2_CLOUD_filename = L2_CLOUD_filename
            self.L1B_GEO_filename = L1B_GEO_filename
            self.L1B_RAD_filename = L1B_RAD_filename
            self.orbit = orbit
            self.scene = scene
            self.production_datetime = production_datetime
            self.build = build
            self.product_counter = product_counter
            self.granule_ID = granule_ID
            self.L2G_granule_ID = L1CG_RAD_granule_ID
            self.L2G_granule_ID = L1CG_RAD_granule_ID
            self.L1CG_RAD_filename = L1CG_RAD_filename
            self.L2G_LSTE_filename = L2G_LSTE_filename
            self.L2G_CLOUD_filename = L2G_CLOUD_filename
            self.PGE_name = PGE_name
            self.PGE_version = PGE_version
        except MissingRunConfigValue as e:
            raise e
        except ECOSTRESSExitCodeException as e:
            raise e
        except Exception as e:
            logger.exception(e)
            raise UnableToParseRunConfig(f"unable to parse run-config file: {filename}")


def L1_L2_RAD_LSTE(
        runconfig_filename: str,
        tiles: List[str] = None,
        variables: List[str] = None,
        process_tiles: bool = True,
        cell_size_degrees: float = CELL_SIZE_DEGREES,
        cell_size_meters: float = CELL_SIZE_METERS,
        gridded_geometry: RasterGrid = None,
        kd_tree: KDTree = None,
        scan_kd_trees: List[KDTree] = None,
        kd_tree_path: str = None,
        search_radius_meters: float = SEARCH_RADIUS_METERS,
        overlap_strategy: str = OVERLAP_STRATEGY,
        projection_system: str = PROJECTION_SYSTEM,
        strip_console: bool = STRIP_CONSOLE) -> int:
    """
    ECOSTRESS Collection 2 L2G L2T LSTE PGE
    :param runconfig_filename: filename for XML run-config
    :param log_filename: filename for logger output
    :return: exit code number
    """
    exit_code = SUCCESS_EXIT_CODE

    try:
        runconfig = L2GL2TRADLSTEConfig(runconfig_filename)
        working_directory = runconfig.working_directory
        granule_ID = runconfig.granule_ID
        log_filename = join(working_directory, "log", f"{granule_ID}.log")
        cl.configure(filename=log_filename, strip_console=strip_console)

        logger.info(f"L1_L2_RAD_LSTE PGE ({cl.val(runconfig.PGE_version)})")
        logger.info(f"L1_L2_RAD_LSTE run-config: {cl.file(runconfig_filename)}")

        logger.info(f"working_directory: {cl.dir(working_directory)}")
        output_directory = runconfig.output_directory
        logger.info(f"output directory: {cl.dir(output_directory)}")
        logger.info(f"log: {cl.file(log_filename)}")
        orbit = runconfig.orbit
        logger.info(f"orbit: {cl.val(orbit)}")
        scene = runconfig.scene
        logger.info(f"scene: {cl.val(scene)}")
        build = runconfig.build
        logger.info(f"build: {cl.val(build)}")
        product_counter = runconfig.product_counter
        logger.info(f"product counter: {cl.val(product_counter)}")
        L2_LSTE_filename = runconfig.L2_LSTE_filename
        logger.info(f"L2_LSTE file: {cl.file(L2_LSTE_filename)}")
        L2_CLOUD_filename = runconfig.L2_CLOUD_filename
        logger.info(f"L2_CLOUD file: {cl.file(L2_CLOUD_filename)}")
        L1B_GEO_filename = runconfig.L1B_GEO_filename
        logger.info(f"L1B_GEO file: {cl.file(L1B_GEO_filename)}")
        L1B_RAD_filename = runconfig.L1B_RAD_filename
        logger.info(f"L1B_RAD file: {cl.file(L1B_RAD_filename)}")

        L1B_RAD_granule = L1BRAD(
            L1B_RAD_filename=L1B_RAD_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L2_LSTE_filename=L2_LSTE_filename,
            L1B_GEO_filename=L1B_GEO_filename
        )

        L2_LSTE_granule = L2LSTE.open(
            L2_LSTE_filename=L2_LSTE_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename
        )

        L2_CLOUD_granule = L2CLOUD.open(
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename
        )

        swath_geometry = L2_LSTE_granule.geometry

        time_UTC = L2_LSTE_granule.time_UTC
        granule_ID = runconfig.granule_ID
        L2G_granule_ID = runconfig.L2G_granule_ID
        L2G_granule_ID = runconfig.L2G_granule_ID
        L1CG_RAD_filename = runconfig.L1CG_RAD_filename
        L2G_LSTE_filename = runconfig.L2G_LSTE_filename
        L2G_CLOUD_filename = runconfig.L2G_CLOUD_filename
        PGE_name = "L1_L2_RAD_LSTE"
        PGE_version = ECOSTRESS.PGEVersion

        if kd_tree_path is None:
            kd_tree_path = join(working_directory, f"{granule_ID}.kdtree")

        input_filenames = [
            L2_LSTE_filename,
            L2_CLOUD_filename,
            L1B_GEO_filename,
            L1B_RAD_filename
        ]

        L1CG_RAD_input_filenames = [
            L2_CLOUD_filename,
            L1B_GEO_filename,
            L1B_RAD_filename
        ]

        L2G_LSTE_input_filenames = [
            L2_LSTE_filename,
            L2_CLOUD_filename,
            L1B_GEO_filename,
            L1B_RAD_filename
        ]

        L2G_CLOUD_input_filenames = [
            L2_CLOUD_filename,
            L1B_GEO_filename,
        ]

        L1CG_RAD_browse_filename = L1CG_RAD_filename.replace(".h5", ".png")

        if gridded_geometry is None:
            if projection_system == "local_UTM":
                gridded_geometry = swath_geometry.UTM(cell_size_meters)
            elif projection_system == "global_geographic":
                gridded_geometry = swath_geometry.geographic(cell_size_degrees)
            else:
                raise ValueError(f"unrecognized projection system: {projection_system}")

        if kd_tree_path is not None and not exists(kd_tree_path):
            logger.warning(f"K-D tree file not found: {kd_tree_path}")

        # TODO consolidate redundancy in preparing K-D trees between L1_RAD.py, L2_LSTE.py, L2G_CLOUD.py, L1_L1_RAD_LSTE.py

        if overlap_strategy == "checkerboard":
            logger.info("using checkerboard overlap strategy")

            # treating K-D tree path as filename for whole-scene processing
            kd_tree_path = kd_tree_path

            if kd_tree_path is not None and exists(kd_tree_path):
                logger.info(f"started loading checkerboard K-D tree: {kd_tree_path}")
                timer = Timer()
                kd_tree = KDTree.load(kd_tree_path)
                logger.info(f"finished loading checkerboard K-D tree: {kd_tree_path} ({timer})")
            else:
                logger.info("started building checkerboard K-D tree")
                timer = Timer()

                kd_tree = KDTree(
                    source_geometry=swath_geometry,
                    target_geometry=gridded_geometry,
                    radius_of_influence=search_radius_meters
                )

                logger.info(f"finished building checkerboard K-D tree ({timer})")

                if kd_tree_path is not None:
                    logger.info(f"started saving checkerboard K-D tree: {kd_tree_path}")
                    timer = Timer()
                    kd_tree.save(kd_tree_path)
                    logger.info(f"finished saving checkerboard K-D tree ({timer}): {kd_tree_path}")

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
                            kd_tree_path = join(kd_tree_directory, f"{i:02d}.kdtree")
                            kd_tree.save(kd_tree_path)

        elif overlap_strategy == "remove_105_128":
            logger.info("using clipped-tails overlap strategy")

            kd_tree_path = kd_tree_path

            logger.info(f"original swath shape: {swath_geometry.shape}")
            swath_geometry = clip_tails(swath_geometry)
            logger.info(f"clipped swath shape: {swath_geometry.shape}")

            if kd_tree_path is not None and exists(kd_tree_path):
                logger.info(f"started loading checkerboard K-D tree: {cl.file(kd_tree_path)}")
                timer = Timer()
                kd_tree = KDTree.load(kd_tree_path)
                logger.info(f"finished loading checkerboard K-D tree: {cl.file(kd_tree_path)} ({cl.time(timer)})")
            else:
                logger.info("started building clipped-tails K-D tree")
                timer = Timer()

                kd_tree = KDTree(
                    source_geometry=swath_geometry,
                    target_geometry=gridded_geometry,
                    radius_of_influence=search_radius_meters
                )

                logger.info(f"finished building clipped-tails K-D tree ({cl.time(timer)})")

                if kd_tree_path is not None:
                    logger.info(f"started saving clipped-tails K-D tree: {cl.file(kd_tree_path)}")
                    timer = Timer()
                    kd_tree.save(kd_tree_path)
                    logger.info(
                        f"finished saving clipped-tails K-D tree ({cl.time(timer)}): {cl.file(kd_tree_path)}")
        else:
            raise ValueError(f"unrecognized overlap strategy: {overlap_strategy}")

        gridded_geometry = swath_geometry.geographic(cell_size_degrees)

        if exists(L1CG_RAD_filename) and exists(L1CG_RAD_browse_filename):
            logger.info(f"found L1CG RAD product file: {cl.file(L1CG_RAD_filename)}")
            # TODO there needs to be a file integrity verification here for the previously generated HDF5 file
            L1CG_RAD_granule = L1CGRAD(L1CG_RAD_filename=L1CG_RAD_filename)
        else:
            logger.info(f"generating L1CG RAD product for orbit {cl.val(orbit)} scene {cl.val(scene)}")
            L1CG_RAD_granule = L1CGRAD.from_swath(
                swath_granule=L1B_RAD_granule,
                output_filename=L1CG_RAD_filename,
                variables=variables,
                PGE_name=PGE_name,
                PGE_version=PGE_version,
                build=build,
                input_filenames=L1CG_RAD_input_filenames,
                gridded_geometry=gridded_geometry,
                cell_size_degrees=cell_size_degrees,
                cell_size_meters=cell_size_meters,
                projection_system=projection_system,
                kd_tree=kd_tree,
                scan_kd_trees=scan_kd_trees,
                kd_tree_path=kd_tree_path
            )

            logger.info(f"generating L2G RAD browse image: {cl.file(L1CG_RAD_browse_filename)}")
            L1CG_RAD_granule.write_browse_image(PNG_filename=L1CG_RAD_browse_filename)

        if process_tiles:
            logger.info(f"generating L1CG RAD tile products for orbit {cl.val(orbit)} scene {cl.val(scene)}")
            L1CG_RAD_granule.to_tiles(
                output_directory=output_directory,
                tiles=tiles,
                variables=variables
            )

        if kd_tree_path is not None and exists(kd_tree_path):
            logger.info(f"re-using K-D tree: {cl.file(kd_tree_path)}")

        L2G_LSTE_browse_filename = L2G_LSTE_filename.replace(".h5", ".png")

        if exists(L2G_LSTE_filename) and exists(L2G_LSTE_browse_filename):
            logger.info(f"found L2G LSTE product file: {cl.file(L2G_LSTE_filename)}")
            logger.info(f"found L2G LSTE browse image: {cl.file(L2G_LSTE_browse_filename)}")
            # TODO there needs to be a file integrity verification here for the previously generated HDF5 file
            L2G_LSTE_granule = L2GLSTE(L2G_LSTE_filename=L2G_LSTE_filename)
        else:
            logger.info(f"generating L2G LSTE gridded product for orbit {cl.val(orbit)} scene {cl.val(scene)}")
            L2G_LSTE_granule = L2GLSTE.from_swath(
                swath_granule=L2_LSTE_granule,
                output_filename=L2G_LSTE_filename,
                variables=variables,
                PGE_name=PGE_name,
                PGE_version=PGE_version,
                build=build,
                input_filenames=L2G_LSTE_input_filenames,
                gridded_geometry=gridded_geometry,
                cell_size_degrees=cell_size_degrees,
                cell_size_meters=cell_size_meters,
                projection_system=projection_system,
                kd_tree=kd_tree,
                scan_kd_trees=scan_kd_trees,
                kd_tree_path=kd_tree_path
            )

            logger.info(f"generating L2G LSTE browse image: {cl.file(L2G_LSTE_browse_filename)}")
            L2G_LSTE_granule.write_browse_image(PNG_filename=L2G_LSTE_browse_filename)

        if process_tiles:
            logger.info(f"generating L2T LSTE tile products for orbit {cl.val(orbit)} scene {cl.val(scene)}")
            L2G_LSTE_granule.to_tiles(
                output_directory=output_directory,
                tiles=tiles,
                variables=variables
            )

        L2G_CLOUD_browse_filename = L2G_CLOUD_filename.replace(".h5", ".png")

        if exists(L2G_CLOUD_filename) and exists(L2G_CLOUD_browse_filename):
            logger.info(f"found L2G CLOUD product file: {cl.file(L2G_CLOUD_filename)}")
            logger.info(f"found L2G CLOUD browse image: {cl.file(L2G_CLOUD_browse_filename)}")
            # TODO there needs to be a file integrity verification here for the previously generated HDF5 file
            L2G_CLOUD_granule = L2GCLOUD(L2G_CLOUD_filename=L2G_CLOUD_filename)
        else:
            logger.info(f"generating L2G CLOUD gridded product for orbit {cl.val(orbit)} scene {cl.val(scene)}")
            L2G_CLOUD_granule = L2GCLOUD.from_swath(
                swath_granule=L2_CLOUD_granule,
                output_filename=L2G_CLOUD_filename,
                variables=variables,
                PGE_name=PGE_name,
                PGE_version=PGE_version,
                build=build,
                input_filenames=L2G_CLOUD_input_filenames,
                gridded_geometry=gridded_geometry,
                cell_size_degrees=cell_size_degrees,
                cell_size_meters=cell_size_meters,
                projection_system=projection_system,
                kd_tree=kd_tree,
                scan_kd_trees=scan_kd_trees,
                kd_tree_path=kd_tree_path
            )

            logger.info(f"generating L2G CLOUD browse image: {cl.file(L2G_CLOUD_browse_filename)}")
            L2G_CLOUD_granule.write_browse_image(PNG_filename=L2G_CLOUD_browse_filename)

        # if kd_tree_path is not None and exists(kd_tree_path):
        #     logger.info(f"removing K-D tree: {cl.file(kd_tree_path)}")
        #     remove(kd_tree_path)

    except ECOSTRESSExitCodeException as exception:
        logger.exception(exception)
        exit_code = exception.exit_code

    return exit_code


def main(argv=sys.argv):
    if len(argv) == 1 or "--version" in argv:
        print(f"L1_L2_RAD_LSTE PGE ({ECOSTRESS.PGEVersion})")
        print(f"usage: L1_L2_RAD_LSTE RunConfig.xml")

        if "--version" in argv:
            return SUCCESS_EXIT_CODE
        else:
            return RUNCONFIG_FILENAME_NOT_SUPPLIED

    strip_console = "--strip-console" in argv
    runconfig_filename = str(argv[1])
    exit_code = L1_L2_RAD_LSTE(runconfig_filename=runconfig_filename, strip_console=strip_console)
    logger.info(f"L1_L2_RAD_LSTE exit code: {exit_code}")

    return exit_code


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
