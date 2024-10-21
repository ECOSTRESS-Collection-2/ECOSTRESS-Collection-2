import logging
from shutil import rmtree
import sys
from dateutil import parser
from os.path import exists
from os.path import splitext
from os.path import basename
from os.path import join
from ECOSTRESS import find_ECOSTRESS_C1_files
from ECOSTRESS.L1B_GEO import L1BGEO
from ECOSTRESS.exit_codes import UncalibratedGeolocation
from L2_LSTE.L2_LSTE import L2_LSTE_from_L1B
import colored_logging as cl

logger = logging.getLogger(__name__)

def main(argv=sys.argv):
    if len(argv) == 1:
        print("ECOv002 --date yyyy-mm-dd --orbit OOOOO --scene SSS [--tiles TTTTT,TTTTT,...] --output directory")

    if "--orbit" in argv:
        orbit = int(argv[argv.index("--orbit") + 1])
    else:
        orbit = None

    if "--scene" in argv:
        scene = int(argv[argv.index("--scene") + 1])
    else:
        scene = None

    if "--date" in argv:
        date_UTC = parser.parse(argv[argv.index("--date") + 1]).date()
    else:
        date_UTC = None

    if "--tiles" in argv:
        tiles = argv[argv.index("--tiles") + 1].split(",")
    else:
        tiles = None

    if "--output" in argv:
        output_directory = str(argv[argv.index("--output") + 1])
    else:
        output_directory = None

    try:
        filenames = find_ECOSTRESS_C1_files(
            orbit=orbit,
            scene=scene,
            date_UTC=date_UTC
        )
    except Exception as e:
        logger.warning(f"orbit {orbit} scene {scene} not found")
        return

    if "L1B_RAD" not in filenames or "L1B_GEO" not in filenames:
        logger.warning(f"orbit {orbit} scene {scene} not found")
        return

    L1B_RAD_filename = filenames["L1B_RAD"]
    L1B_GEO_filename = filenames["L1B_GEO"]
    logger.info(f"L1B RAD file: {cl.file(L1B_RAD_filename)}")
    logger.info(f"L1B GEO file: {cl.file(L1B_GEO_filename)}")
    timestamp = splitext(basename(L1B_GEO_filename))[0].split("_")[-3]
    granule_ID = f"ECOv002_L2G_LSTE_{orbit:05d}_{scene:03d}_{timestamp}_0700_01"
    logger.info(f"granule ID: {cl.val(granule_ID)}")
    granule_run_directory = join(output_directory, granule_ID)
    logger.info(f"granule run directory: {cl.dir(granule_run_directory)}")
    ST_filename = join(output_directory, f"{granule_ID}_ST.tif")
    logger.info(f"output surface temperature file: {cl.file(ST_filename)}")
    water_filename = join(output_directory, f"{granule_ID}_water.tif")
    logger.info(f"output water mask file: {cl.file(water_filename)}")

    if exists(ST_filename) and exists(water_filename):
        logger.info(f"output surface temperature file already exists: {cl.file(ST_filename)}")
        logger.info(f"output water mask file already exists: {cl.file(water_filename)}")
        return

    try:
        L1B_GEO_granule = L1BGEO(L1B_GEO_filename=L1B_GEO_filename)
        geometry = L1B_GEO_granule.geometry
    except UncalibratedGeolocation:
        logger.warning(f"orbit {orbit} scene {scene} is uncalibrated")
        return
    except Exception as e:
        logger.warning(f"unable to load geolocation for orbit {orbit} scene {scene}")
        return

    L2_LSTE_granule = L2_LSTE_from_L1B(
        L1B_RAD_filename=L1B_RAD_filename,
        L1B_GEO_filename=L1B_GEO_filename,
        ST_filename=ST_filename,
        water_filename=water_filename,
        working_directory=granule_run_directory,
        use_docker=True
    )

    logger.info(f"removing granule run directory: {cl.dir(granule_run_directory)}")
    rmtree(granule_run_directory, ignore_errors=True)


if __name__ == "__main__":
    main(argv=sys.argv)
