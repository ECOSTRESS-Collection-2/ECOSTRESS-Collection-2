from glob import glob
from os.path import basename, splitext, dirname
from os.path import join

from dateutil import parser


def parse_ECOSTRESSC1_product(ECOSTRESSC1_filename):
    return "_".join(splitext(basename(ECOSTRESSC1_filename))[0].split("_")[1:-5])


def parse_ECOSTRESSC1_orbit(ECOSTRESSC1_filename):
    return int(splitext(basename(ECOSTRESSC1_filename))[-2].split('_')[-5])


def parse_ECOSTRESSC1_scene_ID(ECOSTRESSC1_filename):
    return int(splitext(basename(ECOSTRESSC1_filename))[-2].split('_')[-4])


def parse_ECOSTRESSC1_build(ECOSTRESSC1_filename):
    return splitext(basename(ECOSTRESSC1_filename))[-2].split('_')[-2]


def parse_ECOSTRESSC1_process_count(ECOSTRESSC1_filename):
    return int(splitext(basename(ECOSTRESSC1_filename))[-2].split('_')[-1])


def parse_ECOSTRESSC1_timestamp(ECOSTRESSC1_filename):
    return parser.parse(splitext(basename(ECOSTRESSC1_filename))[-2].split('_')[-3])


def parse_ECOSTRESSC1_filename(ecostress_filename):
    time = parse_ECOSTRESSC1_timestamp(ecostress_filename)
    scene = parse_ECOSTRESSC1_scene_ID(ecostress_filename)
    orbit = parse_ECOSTRESSC1_orbit(ecostress_filename)
    build = parse_ECOSTRESSC1_build(ecostress_filename)
    process_count = parse_ECOSTRESSC1_process_count(ecostress_filename)

    return orbit, scene, time, build, process_count


def find_corresponding_filename(given_filename: str, product: str) -> str:
    orbit, scene, time, build, process_count = parse_ECOSTRESSC1_filename(given_filename)
    directory = dirname(given_filename)
    filename_pattern = f"ECO*_{product}_{orbit:05d}_{scene:03d}_*.h5"
    filenames = sorted(glob(join(directory, filename_pattern)))

    if len(filenames) == 0:
        raise IOError(f"no matching filename for product {product} from filename {given_filename}")

    filename = filenames[-1]

    return filename
