import logging
from glob import glob
from os.path import basename, splitext, dirname
from os.path import join

from dateutil import parser
import colored_logging as cl

logger = logging.getLogger(__name__)

def parse_ECOSTRESS_product(ECOSTRESS_filename):
    return "_".join(splitext(basename(ECOSTRESS_filename))[0].split("_")[1:-5])


def parse_ECOSTRESS_orbit(ECOSTRESS_filename):
    if ECOSTRESS_filename.endswith(".h5"):
        return int(splitext(basename(ECOSTRESS_filename))[-2].split('_')[-5])
    elif ECOSTRESS_filename.endswith(".zip"):
        return int(splitext(basename(ECOSTRESS_filename))[-2].split('_')[-6])


def parse_ECOSTRESS_scene_ID(ECOSTRESS_filename):
    if ECOSTRESS_filename.endswith(".h5"):
        return int(splitext(basename(ECOSTRESS_filename))[-2].split('_')[-4])
    elif ECOSTRESS_filename.endswith(".zip"):
        return int(splitext(basename(ECOSTRESS_filename))[-2].split('_')[-5])


def parse_ECOSTRESS_build(ECOSTRESS_filename):
    return splitext(basename(ECOSTRESS_filename))[-2].split('_')[-2]


def parse_ECOSTRESS_process_count(ECOSTRESS_filename):
    return int(splitext(basename(ECOSTRESS_filename))[-2].split('_')[-1])


def parse_ECOSTRESS_timestamp(ECOSTRESS_filename):
    return parser.parse(splitext(basename(ECOSTRESS_filename))[-2].split('_')[-3])


def parse_ECOSTRESS_filename(ecostress_filename):
    time = parse_ECOSTRESS_timestamp(ecostress_filename)
    scene = parse_ECOSTRESS_scene_ID(ecostress_filename)
    orbit = parse_ECOSTRESS_orbit(ecostress_filename)
    build = parse_ECOSTRESS_build(ecostress_filename)
    process_count = parse_ECOSTRESS_process_count(ecostress_filename)

    return orbit, scene, time, build, process_count


def find_corresponding_filename(given_filename: str, product: str, match_build: bool = False) -> str:
    orbit, scene, time, build, process_count = parse_ECOSTRESS_filename(given_filename)
    directory = dirname(given_filename)

    if not match_build:
        build = "*"

    filename_pattern = join(directory, f"ECO*_{product}_{orbit:05d}_{scene:03d}_*_{build}_*.h5")
    logger.info(f"searching pattern: {cl.file(filename_pattern)}")
    filenames = sorted(glob(filename_pattern))

    if len(filenames) == 0:
        raise IOError(f"no matching filename for product {product} from filename {given_filename}")

    filename = filenames[-1]
    logger.info(f"found: {cl.file(filename)}")

    return filename
