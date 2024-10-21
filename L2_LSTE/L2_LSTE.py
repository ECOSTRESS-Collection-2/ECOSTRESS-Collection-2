import logging
import posixpath
import shutil
import socket
import sys
from datetime import datetime
from glob import glob
from os import makedirs, system, chmod, rename
from os.path import join, abspath, dirname, expanduser, exists, basename, splitext
from shutil import which
from stat import S_IRWXU, S_IRWXG, S_IRWXO
from uuid import uuid4

from dateutil import parser

import ECOSTRESS
import colored_logging as cl
from ECOSTRESS.L2_LSTE import ECOv002L2LSTE, L2LSTE
from ECOSTRESS.exit_codes import SUCCESS_EXIT_CODE, RUNCONFIG_FILENAME_NOT_SUPPLIED, \
    MissingRunConfigValue, UnableToParseRunConfig, ECOSTRESSExitCodeException
from ECOSTRESS.find_ECOSTRESS_C1_scene import find_ECOSTRESS_C1_scene
from ECOSTRESS.runconfig import ECOSTRESSRunConfig

with open(join(abspath(dirname(__file__)), "version.txt")) as f:
    version = f.read()

__version__ = version

L2_LSTE_TEMPLATE = join(abspath(dirname(__file__)), "L2_LSTE.xml")
DEFAULT_BUILD = "0700"
DEFAULT_WORKING_DIRECTORY = "."
DEFAULT_DOCKER_INPUT_DIRECTORY = "/input"
DEFAULT_NATIVE_NWP_DIRECTORY = "/NWP"
DEFAULT_DOCKER_NWP_DIRECTORY = "/NWP"
DEFAULT_OUTPUT_DIRECTORY = "L2_LSTE_output"
DEFAULT_DOCKER_OUTPUT_DIRECTORY = "/output"
DEFAULT_OSP_DIRECTORY = "/build/ecostress-level2/OSP"
DEFAULT_NWP_DIRECTORY = "/NWP/GEOS5"
DEFAULT_EXECUTABLE_FILENAME = "docker run ecostress/ecostress-level2 /pge/L2_PGE.sh"

logger = logging.getLogger(__name__)


def generate_L2_LSTE_runconfig(
        L1B_RAD_filename: str,
        L1B_GEO_filename: str = None,
        orbit: int = None,
        scene: int = None,
        working_directory: str = None,
        NWP_directory: str = None,
        OSP_directory: str = None,
        executable_filename: str = None,
        native_input_directory: str = None,
        internal_input_directory: str = None,
        native_output_directory: str = None,
        internal_output_directory: str = None,
        runconfig_filename: str = None,
        log_filename: str = None,
        build: str = None,
        processing_node: str = None,
        production_datetime: datetime = None,
        job_ID: str = None,
        instance_ID: str = None,
        product_counter: int = None,
        template_filename: str = None,
        use_docker: bool = False) -> str:
    L1B_RAD_filename = abspath(expanduser(L1B_RAD_filename))

    if not exists(L1B_RAD_filename):
        raise IOError(f"L1B RAD file not found: {L1B_RAD_filename}")

    logger.info(f"L1B RAD file: {cl.file(L1B_RAD_filename)}")
    source_granule_ID = splitext(basename(L1B_RAD_filename))[0]
    logger.info(f"source granule ID: {cl.name(source_granule_ID)}")

    if orbit is None:
        orbit = int(source_granule_ID.split("_")[-5])

    logger.info(f"orbit: {cl.val(orbit)}")

    if scene is None:
        scene = int(source_granule_ID.split("_")[-4])

    logger.info(f"scene: {cl.val(scene)}")

    if L1B_GEO_filename is None:
        directory = abspath(expanduser(dirname(L1B_RAD_filename)))
        pattern = join(directory, f"*_L1B_GEO_{orbit:05d}_{scene:03d}_*.h5")
        logger.info(f"searching for L1B GEO: {cl.val(pattern)}")
        candidates = sorted(glob(pattern))

        if len(candidates) == 0:
            raise ValueError("no L1B GEO filename given or found")

        L1B_GEO_filename = candidates[-1]

    logger.info(f"L1B GEO file: {cl.file(L1B_GEO_filename)}")

    if template_filename is None:
        template_filename = L2_LSTE_TEMPLATE

    template_filename = abspath(expanduser(template_filename))

    if executable_filename is None:
        executable_filename = which("L2_LSTE")

    if executable_filename is None:
        # raise EnvironmentError("executable L2_LSTE not found")
        executable_filename = DEFAULT_EXECUTABLE_FILENAME

    if native_input_directory is None:
        native_input_directory = working_directory

    logger.info(f"native input directory: {native_input_directory}")

    if internal_input_directory is None:
        if use_docker:
            internal_input_directory = DEFAULT_DOCKER_INPUT_DIRECTORY
        else:
            internal_input_directory = native_input_directory

    logger.info(f"internal input directory: {internal_input_directory}")

    if native_output_directory is None:
        native_output_directory = join(working_directory, DEFAULT_OUTPUT_DIRECTORY)

    native_output_directory = abspath(expanduser(native_output_directory))

    logger.info(f"native output directory: {native_output_directory}")

    if internal_output_directory is None:
        if use_docker:
            internal_output_directory = join(working_directory, DEFAULT_OUTPUT_DIRECTORY)
        else:
            internal_output_directory = native_output_directory

    internal_output_directory = abspath(expanduser(internal_output_directory))

    logger.info(f"internal output directory: {internal_output_directory}")

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

    L1B_GEO_filename = abspath(expanduser(L1B_GEO_filename))
    timestamp = splitext(basename(L1B_GEO_filename))[0].split("_")[-3]
    granule_ID = f"ECOv002_L2_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"

    if runconfig_filename is None:
        runconfig_filename = join(native_output_directory, f"{granule_ID}.xml")

    runconfig_filename = abspath(expanduser(runconfig_filename))

    logger.info(f"run-config file: {runconfig_filename}")

    if log_filename is None:
        log_filename = join(working_directory, "log", f"{granule_ID}.log")

    log_filename = abspath(expanduser(log_filename))

    if working_directory is None:
        working_directory = granule_ID

    working_directory = abspath(expanduser(working_directory))

    if NWP_directory is None:
        NWP_directory = DEFAULT_NWP_DIRECTORY

    if OSP_directory is None:
        OSP_directory = DEFAULT_OSP_DIRECTORY

    logger.info(f"generating run-config for orbit {cl.val(orbit)} scene {cl.val(scene)}")
    logger.info(f"loading L2_LSTE template: {cl.file(template_filename)}")

    with open(template_filename, "r") as file:
        template = file.read()

    logger.info(f"orbit: {cl.val(orbit)}")
    template = template.replace("orbit_number", f"{orbit:05d}")
    logger.info(f"scene: {cl.val(scene)}")
    template = template.replace("scene_ID", f"{scene:03d}")

    if use_docker:
        internal_L1B_GEO_filename = join(internal_input_directory, basename(L1B_GEO_filename))
    else:
        internal_L1B_GEO_filename = L1B_GEO_filename

    logger.info(f"L1B_GEO file: {cl.file(internal_L1B_GEO_filename)}")
    template = template.replace("L1B_GEO_filename", internal_L1B_GEO_filename)

    if use_docker:
        internal_L1B_RAD_filename = join(internal_input_directory, basename(L1B_RAD_filename))
    else:
        internal_L1B_RAD_filename = L1B_RAD_filename

    logger.info(f"L1B_RAD file: {cl.file(internal_L1B_RAD_filename)}")
    template = template.replace("L1B_RAD_filename", internal_L1B_RAD_filename)
    logger.info(f"NWP directory: {cl.dir(NWP_directory)}")
    template = template.replace("NWP_directory", NWP_directory)
    logger.info(f"OSP directory: {cl.dir(OSP_directory)}")
    template = template.replace("OSP_directory", OSP_directory)
    logger.info(f"executable: {cl.file(executable_filename)}")
    template = template.replace("executable_filename", executable_filename)
    logger.info(f"output directory: {cl.dir(internal_output_directory)}")
    template = template.replace("output_directory", internal_output_directory)
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


def L2_LSTE_runconfig_from_C1(orbit: int, scene: int, runconfig_filename: str, use_docker: bool = False) -> str:
    filenames = find_ECOSTRESS_C1_scene(
        orbit=orbit,
        scene=scene
    )

    L2_LSTE_filename = filenames["L2_LSTE"]
    L2_CLOUD_filename = filenames["L2_CLOUD"]
    L1B_GEO_filename = filenames["L1B_GEO"]

    runconfig = generate_L2_LSTE_runconfig(
        orbit=orbit,
        scene=scene,
        L2_LSTE_filename=L2_LSTE_filename,
        L2_CLOUD_filename=L2_CLOUD_filename,
        L1B_GEO_filename=L1B_GEO_filename,
        runconfig_filename=runconfig_filename,
        use_docker=use_docker
    )

    return runconfig


class L2GL2TRADLSTEConfig(ECOSTRESSRunConfig):
    def __init__(self, filename: str):
        logger.info(f"loading L2_LSTE run-config: {cl.file(filename)}")
        runconfig = self.read_runconfig(filename)

        try:
            if "StaticAncillaryFileGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing StaticAncillaryFileGroup in L2_LSTE run-config: {filename}")

            if "ProductPathGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup in L2_LSTE run-config: {filename}")

            if "ProductPath" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductPath in L2_LSTE run-config: {filename}")

            output_directory = abspath(runconfig["ProductPathGroup"]["ProductPath"])
            logger.info(f"output directory: {cl.dir(output_directory)}")

            if "InputFileGroup" not in runconfig:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup in L2_LSTE run-config: {filename}")

            if "L1B_GEO" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L1B_GEO in L2_LSTE run-config: {filename}")

            L1B_GEO_filename = abspath(runconfig["InputFileGroup"]["L1B_GEO"])
            logger.info(f"L1B_GEO file: {cl.file(L1B_GEO_filename)}")

            if "L1B_RAD" not in runconfig["InputFileGroup"]:
                raise MissingRunConfigValue(
                    f"missing InputFileGroup/L1B_RAD in L2_LSTE run-config: {filename}")

            L1B_RAD_filename = abspath(runconfig["InputFileGroup"]["L1B_RAD"])
            logger.info(f"L1B_RAD file: {cl.file(L1B_RAD_filename)}")

            orbit = int(runconfig["Geometry"]["OrbitNumber"])
            logger.info(f"orbit: {cl.val(orbit)}")

            if "SceneId" not in runconfig["Geometry"]:
                raise MissingRunConfigValue(
                    f"missing Geometry/SceneId in L2_LSTE run-config: {filename}")

            scene = int(runconfig["Geometry"]["SceneId"])
            logger.info(f"scene: {cl.val(scene)}")

            if "ProductionDateTime" not in runconfig["JobIdentification"]:
                raise MissingRunConfigValue(
                    f"missing JobIdentification/ProductionDateTime in L2_LSTE run-config {filename}")

            production_datetime = parser.parse(runconfig["JobIdentification"]["ProductionDateTime"])
            logger.info(f"production time: {cl.time(production_datetime)}")

            if "BuildID" not in runconfig["PrimaryExecutable"]:
                raise MissingRunConfigValue(
                    f"missing PrimaryExecutable/BuildID in L2_LSTE run-config {filename}")

            build = str(runconfig["PrimaryExecutable"]["BuildID"])

            if "FullPathname" not in runconfig["PrimaryExecutable"]:
                raise MissingRunConfigValue(
                    f"missing PrimaryExecutable/FullPathname in L2_LSTE run-config {filename}")

            executable_filename = str(runconfig["PrimaryExecutable"]["FullPathname"])

            if "ProductCounter" not in runconfig["ProductPathGroup"]:
                raise MissingRunConfigValue(
                    f"missing ProductPathGroup/ProductCounter in L2_LSTE run-config {filename}")

            PGE_name = "L2_LSTE"
            PGE_version = ECOSTRESS.PGEVersion
            product_counter = int(runconfig["ProductPathGroup"]["ProductCounter"])
            time_UTC = datetime.strptime(basename(L1B_GEO_filename).split("_")[-3], "%Y%m%dT%H%M%S")
            timestamp = f"{time_UTC:%Y%m%dT%H%M%S}"

            granule_ID = f"ECOv002_L2_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L2_LSTE_granule_ID = f"ECOv002_L2_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L2_CLOUD_granule_ID = f"ECOv002_L2_CLOUD_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L1CG_RAD_granule_ID = f"ECOv002_L1CG_RAD_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"
            L2G_LSTE_granule_ID = f"ECOv002_L2G_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_{build}_{product_counter:02d}"

            L2_LSTE_filename = join(output_directory, f"{L2_LSTE_granule_ID}.h5")
            L2_CLOUD_filename = join(output_directory, f"{L2_CLOUD_granule_ID}.h5")
            L1CG_RAD_filename = join(output_directory, f"{L1CG_RAD_granule_ID}.h5")
            L2G_LSTE_filename = join(output_directory, f"{L2G_LSTE_granule_ID}.h5")

            self.output_directory = output_directory
            self.orbit = orbit
            self.scene = scene
            self.production_datetime = production_datetime
            self.build = build
            self.executable_filename = executable_filename
            self.product_counter = product_counter
            self.granule_ID = granule_ID

            self.L2_LSTE_granule_ID = L2_LSTE_granule_ID
            self.L2_CLOUD_granule_ID = L2_CLOUD_granule_ID
            self.L1CG_RAD_granule_ID = L1CG_RAD_granule_ID
            self.L2G_LSTE_granule_ID = L2G_LSTE_granule_ID

            self.L1B_GEO_filename = L1B_GEO_filename
            self.L1B_RAD_filename = L1B_RAD_filename
            self.L2_LSTE_filename = L2_LSTE_filename
            self.L2_CLOUD_filename = L2_CLOUD_filename
            self.L1CG_RAD_filename = L1CG_RAD_filename
            self.L2G_LSTE_filename = L2G_LSTE_filename

            self.PGE_name = PGE_name
            self.PGE_version = PGE_version
        except MissingRunConfigValue as e:
            raise e
        except ECOSTRESSExitCodeException as e:
            raise e
        except Exception as e:
            logger.exception(e)
            raise UnableToParseRunConfig(f"unable to parse run-config file: {filename}")


def L2_LSTE(
        runconfig_filename: str,
        working_directory: str = None,
        native_input_directory: str = None,
        internal_input_directory: str = None,
        native_NWP_directory: str = None,
        internal_NWP_directory: str = None,
        native_output_directory: str = None,
        internal_output_directory: str = None,
        executable_filename: str = None,
        use_docker: bool = False) -> L2LSTE:
    """
    ECOSTRESS Collection 2 L2G L2T LSTE PGE
    :param runconfig_filename: filename for XML run-config
    :param log_filename: filename for logger output
    :return: exit code number
    """
    if working_directory is None:
        working_directory = DEFAULT_WORKING_DIRECTORY

    if native_input_directory is None:
        native_input_directory = working_directory

    if internal_input_directory is None:
        if use_docker:
            internal_input_directory = DEFAULT_DOCKER_INPUT_DIRECTORY
        else:
            internal_input_directory = native_input_directory

    if native_NWP_directory is None:
        native_NWP_directory = DEFAULT_NATIVE_NWP_DIRECTORY

    if internal_NWP_directory is None:
        if use_docker:
            internal_NWP_directory = DEFAULT_DOCKER_NWP_DIRECTORY
        else:
            internal_NWP_directory = native_NWP_directory

    if native_output_directory is None:
        native_output_directory = join(working_directory, DEFAULT_OUTPUT_DIRECTORY)

    if internal_output_directory is None:
        if use_docker:
            internal_output_directory = DEFAULT_DOCKER_OUTPUT_DIRECTORY
        else:
            internal_output_directory = native_output_directory

    if executable_filename is None:
        executable_filename = join(abspath(dirname(__file__)), "ecostress-level2", "src", "L2_PGE")

    logger.info(f"L2_LSTE PGE ({cl.val(__version__)})")
    logger.info(f"L2_LSTE run-config: {cl.file(runconfig_filename)}")
    internal_runconfig_filename = posixpath.join(internal_output_directory, basename(runconfig_filename))
    native_input_directory = abspath(expanduser(native_input_directory))
    native_output_directory = abspath(expanduser(native_output_directory))
    chmod(native_output_directory, S_IRWXU | S_IRWXG | S_IRWXO)
    native_NWP_directory = abspath(expanduser(native_NWP_directory))
    working_directory = abspath(expanduser(working_directory))
    makedirs(native_output_directory, exist_ok=True)

    runconfig = L2GL2TRADLSTEConfig(runconfig_filename)
    L1B_GEO_filename = join(native_input_directory, basename(runconfig.L1B_GEO_filename))
    L2_LSTE_filename = join(native_output_directory, basename(runconfig.L2_LSTE_filename))
    L2_CLOUD_filename = join(native_output_directory, basename(runconfig.L2_CLOUD_filename))

    if exists(L2_LSTE_filename):
        logger.info(f"L2 LSTE already processed: {cl.file(L2_LSTE_filename)}")
    else:
        logger.info(f"L2 LSTE not found: {cl.file(L2_LSTE_filename)}")

    if exists(L2_CLOUD_filename):
        logger.info(f"L2 CLOUD already processed: {cl.file(L2_CLOUD_filename)}")
    else:
        logger.info(f"L2 CLOUD not found: {cl.file(L2_CLOUD_filename)}")

    if exists(L2_LSTE_filename) and exists(L2_CLOUD_filename):
        granule = L2LSTE.open(
            L2_LSTE_filename=L2_LSTE_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename
        )

        return granule

    if use_docker:
        logger.info("calling L2 PGE using docker")
        logger.info(f"internal run-config filename: {internal_runconfig_filename}")

        command = "docker run --rm " \
                  f"-v {native_input_directory}:{internal_input_directory} " \
                  f"-v {native_output_directory}:{internal_output_directory} " \
                  f"-v {native_NWP_directory}:{internal_NWP_directory} " \
                  f"-w {working_directory} " \
                  f"ecostress/ecostress-level2 /pge/L2_PGE.sh {internal_runconfig_filename}"
    else:
        logger.info("calling l2 PGE without docker")
        logger.info(f"executable filename: {executable_filename}")
        logger.info(f"run-config filename: {runconfig_filename}")

        command = f"{executable_filename} {runconfig_filename}"

    logger.info(command)
    system(command)

    expanded_L2_LSTE_filename = join(dirname(L2_LSTE_filename), basename(L2_LSTE_filename).replace("ECOv002_", "ECOSTRESS_"))
    
    if exists(expanded_L2_LSTE_filename):
        logger.warning(f"renaming {expanded_L2_LSTE_filename} -> {L2_LSTE_filename}")
        rename(expanded_L2_LSTE_filename, L2_LSTE_filename)

    expanded_L2_CLOUD_filename = join(dirname(L2_CLOUD_filename), basename(L2_CLOUD_filename).replace("ECOv002_", "ECOSTRESS_"))
    
    if exists(expanded_L2_CLOUD_filename):
        logger.warning(f"renaming {expanded_L2_CLOUD_filename} -> {L2_CLOUD_filename}")
        rename(expanded_L2_CLOUD_filename, L2_CLOUD_filename)

    if not exists(L2_LSTE_filename):
        raise IOError(f"L2 LSTE file not found: {L2_LSTE_filename}")

    if not exists(L2_CLOUD_filename):
        raise IOError(f"L2 CLOUD file not found: {L2_CLOUD_filename}")

    logger.info(f"ingested L1B GEO file: {cl.file(L1B_GEO_filename)}")
    logger.info(f"generated L2 LSTE file: {cl.file(L2_LSTE_filename)}")
    logger.info(f"generated L2 CLOUD file: {cl.file(L2_CLOUD_filename)}")

    granule = L2LSTE.open(
        L2_LSTE_filename=L2_LSTE_filename,
        L2_CLOUD_filename=L2_CLOUD_filename,
        L1B_GEO_filename=L1B_GEO_filename
    )

    return granule


def L2_LSTE_from_L1B(
        L1B_RAD_filename: str,
        L1B_GEO_filename: str = None,
        ST_filename: str = None,
        water_filename: str = None,
        working_directory: str = None,
        native_input_directory: str = None,
        internal_input_directory: str = None,
        native_output_directory: str = None,
        internal_output_directory: str = None,
        OSP_directory: str = None,
        use_docker: bool = True) -> L2LSTE:
    logger.info("processing L2 LSTE from L1B")

    if working_directory is None:
        working_directory = dirname(L1B_RAD_filename)

    logger.info(f"working directory: {working_directory}")

    if native_input_directory is None:
        native_input_directory = working_directory

    logger.info(f"native input directory: {native_input_directory}")

    if internal_input_directory is None:
        if use_docker:
            internal_input_directory = DEFAULT_DOCKER_INPUT_DIRECTORY
        else:
            internal_input_directory = native_input_directory

    logger.info(f"internal input directory: {internal_input_directory}")

    if native_output_directory is None:
        native_output_directory = join(working_directory, DEFAULT_OUTPUT_DIRECTORY)

    logger.info(f"native output directory: {native_output_directory}")

    if internal_output_directory is None:
        if use_docker:
            internal_output_directory = DEFAULT_DOCKER_OUTPUT_DIRECTORY
        else:
            internal_output_directory = native_output_directory

    logger.info(f"internal output directory: {internal_output_directory}")

    if OSP_directory is None:
        if use_docker:
            OSP_directory = DEFAULT_OSP_DIRECTORY
        else:
            OSP_directory = join(abspath(dirname(__file__)), "ecostress-level2", "OSP")

    logger.info(f"OSP directory: {OSP_directory}")

    runconfig_filename = generate_L2_LSTE_runconfig(
        L1B_RAD_filename=L1B_RAD_filename,
        L1B_GEO_filename=L1B_GEO_filename,
        working_directory=working_directory,
        native_input_directory=native_input_directory,
        internal_input_directory=internal_input_directory,
        native_output_directory=native_output_directory,
        internal_output_directory=internal_output_directory,
        OSP_directory=OSP_directory,
        use_docker=use_docker
    )

    if use_docker:
        makedirs(native_input_directory, exist_ok=True)

        L1B_GEO_filename_link = join(native_input_directory, basename(L1B_GEO_filename))

        if not exists(L1B_GEO_filename_link):
            shutil.copy(L1B_GEO_filename, L1B_GEO_filename_link)

        L1B_RAD_filename_link = join(native_input_directory, basename(L1B_RAD_filename))

        if not exists(L1B_RAD_filename_link):
            shutil.copy(L1B_RAD_filename, L1B_RAD_filename_link)

    granule = L2_LSTE(
        runconfig_filename=runconfig_filename,
        working_directory=working_directory,
        native_input_directory=native_input_directory,
        internal_input_directory=internal_input_directory,
        native_output_directory=native_output_directory,
        internal_output_directory=internal_output_directory,
        use_docker=use_docker
    )

    if ST_filename is not None:
        if water_filename is None:
            water_filename = ST_filename.replace("_ST.tif", "_water.tif")

        geometry = granule.geometry.geographic(0.0006)
        ST_C = granule.ST_C.to_geometry(geometry)
        water = granule.water.to_geometry(geometry)
        logger.info(f"writing surface temperature to GeoTIFF: {ST_filename}")
        ST_C.to_geotiff(ST_filename)
        logger.info(f"writing water mask to GeoTIFF: {water_filename}")
        water.to_geotiff(water_filename)

    return granule


def main(argv=sys.argv):
    if len(argv) == 1 or "--version" in argv:
        print(f"L2_LSTE PGE ({__version__})")
        print(f"usage: L2_LSTE RunConfig.xml")

        if "--version" in argv:
            return SUCCESS_EXIT_CODE
        else:
            return RUNCONFIG_FILENAME_NOT_SUPPLIED

    runconfig_filename = str(argv[1])
    exit_code = L2_LSTE(runconfig_filename=runconfig_filename)
    logger.info(f"L2_LSTE exit code: {exit_code}")

    return exit_code


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
