import logging
import sys
from datetime import datetime
from glob import glob
from os.path import join
from time import perf_counter
from typing import Union

from dateutil import parser

import colored_logging as cl
from .ECOSTRESS_C1_filenames import parse_ECOSTRESSC1_build
from .ECOSTRESS_C1_filenames import parse_ECOSTRESSC1_product

DEFAULT_ECOSTRESS_STORE_ROOT = "/ops/store*/PRODUCTS"
DEFAULT_BUILD = "0601"


def find_ECOSTRESS_C1_scene(
        orbit: int,
        scene: int,
        date_UTC: Union[datetime, str] = None,
        store_root: str = None,
        build: str = None) -> dict:
    """
    Search for ECOSTRESS files associated with a scene and orbit.
    :param orbit: orbit number
    :param scene: scene number
    :return: dict organized by build then product
    """
    logger = logging.getLogger(__name__)

    if store_root is None:
        store_root = DEFAULT_ECOSTRESS_STORE_ROOT

    if build is None:
        build = DEFAULT_BUILD

    if isinstance(date_UTC, str):
        date_UTC = parser.parse(date_UTC).date()

    orbit = int(orbit)
    scene = int(scene)
    logger.info(f"searching for ECOSTRESS Collection 1 orbit {cl.val(orbit)} scene {cl.val(scene)}")
    file_pattern = f"ECOSTRESS_*_{orbit:05d}_{scene:03d}_*_*_*.h5"

    if date_UTC is not None:
        logger.info(f"searching with date: {cl.time(date_UTC)}")
        year = date_UTC.year
        month = date_UTC.month
        day = date_UTC.day

        path_pattern = join(
            store_root,
            "*",
            f"{year:04d}",
            f"{month:02d}",
            f"{day:02d}",
            file_pattern
        )

    else:
        path_pattern = join(store_root, "**", file_pattern)

    logger.info(f"searching pattern: {cl.val(path_pattern)}")
    start_time = perf_counter()
    filenames = sorted(glob(path_pattern, recursive=True))
    end_time = perf_counter()
    duration = end_time - start_time
    logger.info(f"found {cl.val(len(filenames))} files ({cl.time(f'{duration:0.2f}s')})")

    if len(filenames) == 0:
        raise FileNotFoundError(f"orbit {orbit} scene {scene} not found")

    results = {}

    for filename in filenames:
        file_build = parse_ECOSTRESSC1_build(filename)
        product = parse_ECOSTRESSC1_product(filename)

        if file_build == build:
            results[product] = filename

    # for product in sorted(results.keys()):
    #     logger.info(f"{product}: {results[product]}")

    return results


def main(argv=sys.argv):
    if len(argv) == 1:
        print("find-ECOSTRESS-C1 OOOOO SSS [--date yyyy-mm-dd]")

    cl.configure()
    logger = logging.getLogger(__name__)

    orbit = int(argv[1])
    scene = int(argv[2])

    if "--date" in argv:
        acquisition_date = parser.parse(argv[argv.index("--date") + 1]).date()
    else:
        acquisition_date = None

    find_ECOSTRESS_C1_scene(
        orbit=orbit,
        scene=scene,
        date_UTC=acquisition_date
    )


if __name__ == "__main__":
    main(argv=sys.argv)
