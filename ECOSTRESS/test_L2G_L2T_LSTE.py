import logging
import os
import sys
from datetime import date
from os.path import join, exists

from ECOSTRESS.find_ECOSTRESS_C1_scene import find_ECOSTRESS_C1_scene

import colored_logging as cl
from L1_L2_RAD_LSTE.L1_L2_RAD_LSTE import L1_L2_RAD_LSTE, generate_L1_L2_RAD_LSTE_runconfig

logger = logging.getLogger(__name__)


def test_L1_L2_RAD_LSTE(
        orbit: int,
        scene: int,
        acquisition_date: date = None,
        build: str = None,
        working_directory: str = None,
        output_directory: str = None,
        runconfig_filename: str = None,
        log_filename: str = None):
    if not isinstance(orbit, int):
        orbit = int(orbit)

    if not isinstance(scene, int):
        scene = int(scene)

    if working_directory is None:
        working_directory = join(os.getcwd(), f"L2G_L2T_LSTE_{orbit:05d}_{scene:03d}")

    if output_directory is None:
        output_directory = join(working_directory, "output")

    if runconfig_filename is None:
        runconfig_filename = join(working_directory, f"L2G_L2T_LSTE_{orbit:05d}_{scene:03d}.xml")

    if not exists(runconfig_filename):
        filenames = find_ECOSTRESS_C1_scene(orbit, scene, date_UTC=acquisition_date, build=build)

        if "L2_LSTE" not in filenames:
            raise IOError("L2_LSTE file not found")

        L2_LSTE_filename = filenames["L2_LSTE"]
        logger.info(f"L2_LSTE file: {cl.file(L2_LSTE_filename)}")

        if "L2_CLOUD" not in filenames:
            raise IOError("L2_CLOUD file not found")

        L2_CLOUD_filename = filenames["L2_CLOUD"]
        logger.info(f"L2_CLOUD file: {cl.file(L2_CLOUD_filename)}")

        if "L1B_GEO" not in filenames:
            raise IOError("L1B_GEO file not found")

        L1B_GEO_filename = filenames["L1B_GEO"]
        logger.info(f"L1B_GEO file: {cl.file(L1B_GEO_filename)}")

        if "L1B_RAD" not in filenames:
            raise IOError("L1B_RAD file not found")

        L1B_RAD_filename = filenames["L1B_RAD"]
        logger.info(f"L1B_RAD file: {cl.file(L1B_RAD_filename)}")

        generate_L1_L2_RAD_LSTE_runconfig(
            orbit=orbit,
            scene=scene,
            L2_LSTE_filename=L2_LSTE_filename,
            L2_CLOUD_filename=L2_CLOUD_filename,
            L1B_GEO_filename=L1B_GEO_filename,
            L1B_RAD_filename=L1B_RAD_filename,
            output_directory=output_directory,
            runconfig_filename=runconfig_filename
        )

    L1_L2_RAD_LSTE(runconfig_filename=runconfig_filename, log_filename=log_filename)


def main(argv=sys.argv):
    if len(argv) < 3:
        print("test_L1_L2_RAD_LSTE OOOOO SSS [--date yyyy-mm-dd] [--build BBbb]")
        sys.exit(1)

    orbit = int(argv[1])
    scene = int(argv[2])

    if '--date' in argv:
        acquisition_date = argv[argv.index('--date') + 1]
    else:
        acquisition_date = None

    if '--build' in argv:
        build = argv[argv.index('--build') + 1]
    else:
        build = None

    test_L1_L2_RAD_LSTE(orbit, scene, acquisition_date=acquisition_date, build=build)


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
