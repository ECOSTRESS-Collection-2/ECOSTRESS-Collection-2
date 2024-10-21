from typing import List
import logging
import sys
from datetime import date
from glob import glob
from os import makedirs, symlink
from os.path import join, abspath, expanduser, exists, basename
from shutil import copyfile
from typing import Union

from dateutil import parser

import pandas as pd

import colored_logging as cl
from ECOSTRESS import L2TLSTE
from ECOSTRESS import open_granule
from ECOSTRESS.L1B_GEO import L1BGEO
from ECOSTRESS.exit_codes import ECOSTRESSExitCodeException
from ECOSTRESS.find_ECOSTRESS_C1_geometry import find_ECOSTRESS_C1_geometry
from L1_L2_RAD_LSTE import generate_L1_L2_RAD_LSTE_runconfig, L1_L2_RAD_LSTE, L2GL2TRADLSTEConfig
from L2T_STARS import generate_L2T_STARS_runconfig, L2TSTARSConfig, L2T_STARS
from L2_LSTE import L2_LSTE_from_L1B
from L3T_L4T_ALEXI import generate_L3T_L4T_ALEXI_runconfig, L3TL4TALEXIConfig, L3T_L4T_ALEXI, DEFAULT_ALEXI_SOURCES_DIRECTORY
from L3T_L4T_JET import generate_L3T_L4T_JET_runconfig, L3TL4TJETConfig, L3T_L4T_JET
from sentinel_tile_grid import sentinel_tile_grid

logger = logging.getLogger(__name__)

CENTRAL_STATIC_DIRECTORY = "/home/halverso/sandbox/ECOv002_L3T_L4T/L3T_L4T_STATIC"
CENTRAL_L2T_STARS_SOURCES = "/home/halverso/sandbox/ECOv002_L3T_L4T/L2T_STARS_SOURCES"
CENTRAL_L3T_L4T_JET_SOURCES = "/home/halverso/sandbox/ECOv002_L3T_L4T/L3T_L4T_JET_SOURCES"
CENTRAL_L3T_L4T_ALEXI_SOURCES = "/home/halverso/sandbox/ECOv002_L3T_L4T/L3T_L4T_ALEXI_SOURCES"

def find_calval_product_files(output_directory: str, product_name: str, PGE_name: str, tile: str, date_UTC: Union[date, str] = None) -> List[str]:
    if date_UTC is None:
        date_UTC = "**"
    else:
        date_UTC = parser.parse(str(date_UTC)).strftime("%Y-%m-%d")

    pattern = join(output_directory, PGE_name, date_UTC, f"*{product_name}*_{tile}_*.zip")
    logger.info(f"searching cal/val output with pattern: {pattern}")
    filenames = glob(pattern, recursive=True)
    logger.info(f"found {len(filenames)} files")

    return filenames

def find_JET_inputs(output_directory: str, tile: str) -> pd.DataFrame:
    L2T_LSTE_filenames = find_calval_product_files(
        output_directory=output_directory,
        product_name="L2T_LSTE",
        PGE_name="L1_L2_RAD_LSTE",
        tile=tile
    )

    L2T_STARS_filenames = find_calval_product_files(
        output_directory=output_directory,
        product_name="L2T_STARS",
        PGE_name="L2T_STARS",
        tile=tile
    )

    L2T_LSTE_rows = []

    for filename in L2T_LSTE_filenames:
        granule = open_granule(filename)
        date_UTC = granule.date_UTC
        L2T_LSTE_rows.append([date_UTC, filename])
    
    L2T_LSTE_df = pd.DataFrame(L2T_LSTE_rows, columns=["date_UTC", "L2T_LSTE_filename"])

    L2T_STARS_rows = []

    for filename in L2T_STARS_filenames:
        granule = open_granule(filename)
        date_UTC = granule.date_UTC
        L2T_STARS_rows.append([date_UTC, filename])
    
    L2T_STARS_df = pd.DataFrame(L2T_STARS_rows, columns=["date_UTC", "L2T_STARS_filename"])

    df = L2T_LSTE_df.merge(L2T_STARS_df, on="date_UTC", how="left")
    df.dropna(inplace=True)
    logger.info(f"found {len(df)} matching pairs of L2T LSTE and L2T STARS")

    return df


def reprocess_JET(
        tile: str,
        start_date_UTC: Union[date, str] = None,
        end_date_UTC: Union[date, str] = None,
        start_date_solar: Union[date, str] = None,
        end_date_solar: Union[date, str] = None,
        hour_min: float = None,
        hour_max: float = None,
        max_cloud_percent: float = None,
        results_filename: str = None,
        bboxes_filenames: str = None,
        main_directory: str = None,
        static_directory: str = None,
        SRTM_directory: str = None,
        L2T_STARS_sources_directory: str = None,
        L3T_L4T_JET_sources_directory: str = None,
        ALEXI_directory: str = None,
        output_directory: str = None,
        SZA_cutoff: float = None,
        build: str = "0700",
        halt_with_unhandled_exceptions: bool = False):
    geometry = sentinel_tile_grid.grid(tile)

    if main_directory is None:
        main_directory = abspath(".")

    if SZA_cutoff is None:
        SZA_cutoff = 70

    if static_directory is None:
        static_directory = join(main_directory, "L3T_L4T_STATIC")

    if not exists(static_directory):
        logger.info(f"linking {static_directory} -> {CENTRAL_STATIC_DIRECTORY}")
        symlink(CENTRAL_STATIC_DIRECTORY, static_directory)

    if SRTM_directory is None:
        SRTM_directory = join(main_directory, "SRTM")

    if L2T_STARS_sources_directory is None:
        L2T_STARS_sources_directory = join(main_directory, "L2T_STARS_SOURCES")

    if not exists(static_directory):
        logger.info(f"linking {L2T_STARS_sources_directory} -> {CENTRAL_L2T_STARS_SOURCES}")
        symlink(CENTRAL_L2T_STARS_SOURCES, L2T_STARS_sources_directory)

    if L3T_L4T_JET_sources_directory is None:
        L3T_L4T_JET_sources_directory = join(main_directory, "L3T_L4T_JET_SOURCES")

    if not exists(static_directory):
        logger.info(f"linking {L3T_L4T_JET_sources_directory} -> {CENTRAL_L3T_L4T_JET_SOURCES}")
        symlink(CENTRAL_L3T_L4T_JET_SOURCES, L3T_L4T_JET_sources_directory)

    if ALEXI_directory is None:
        # ALEXI_directory = join(main_directory, "L3T_L4T_ALEXI_SOURCES")
        ALEXI_directory = DEFAULT_ALEXI_SOURCES_DIRECTORY

    if not exists(static_directory):
        logger.info(f"linking {ALEXI_directory} -> {CENTRAL_L3T_L4T_ALEXI_SOURCES}")
        symlink(CENTRAL_L3T_L4T_ALEXI_SOURCES, ALEXI_directory)

    if output_directory is None:
        output_directory = join(main_directory, "output")
    
    JET_inputs_df = find_JET_inputs(
        output_directory=output_directory,
        tile=tile
    )

    for i, (date_UTC, L2T_LSTE_filename, L2T_STARS_filename) in JET_inputs_df.iterrows():
        L2T_LSTE_granule = open_granule(L2T_LSTE_filename)
        orbit = L2T_LSTE_granule.orbit
        scene = L2T_LSTE_granule.scene
        tile = L2T_LSTE_granule.tile
        L3T_L4T_JET_output = join(output_directory, "L3T_L4T_JET", date_UTC.strftime("%Y-%m-%d"))
        time_UTC = L2T_LSTE_granule.time_UTC
        timestamp = f"{time_UTC:%Y%m%dT%H%M%S}"
        run_ID = f"ECOv002_{orbit:05d}_{scene:03d}_{timestamp}"
        run_directory = join(main_directory, "runs", run_ID)

        L3T_JET_filenames = find_calval_product_files(
            output_directory=output_directory,
            product_name="L3T_JET",
            PGE_name="L3T_L4T_JET",
            tile=tile
        )

        if len(L3T_JET_filenames) > 0:
            logger.info(f"found previously completed L3T JET for orbit {orbit} scene {scene} tile {tile} on date {date_UTC}")
            continue

        logger.info(f"processing L3T JET orbit {orbit} scene {scene} tile {tile} on date {date_UTC}")

        try:
            L3T_L4T_JET_runconfig_filename = generate_L3T_L4T_JET_runconfig(
                orbit=orbit,
                scene=scene,
                tile=tile,
                L2T_LSTE_filename=L2T_LSTE_filename,
                L2T_STARS_filename=L2T_STARS_filename,
                working_directory=run_directory,
                sources_directory=L3T_L4T_JET_sources_directory,
                static_directory=static_directory,
                output_directory=L3T_L4T_JET_output,
                SRTM_directory=SRTM_directory
            )

            L3T_L4T_JET_config = L3TL4TJETConfig(L3T_L4T_JET_runconfig_filename)
            L3T_L4T_JET_sources_directory = L3T_L4T_JET_config.output_directory
            logger.info(f"L3T L4T directory: {L3T_L4T_JET_sources_directory}")
            exit_code = L3T_L4T_JET(runconfig_filename=L3T_L4T_JET_runconfig_filename)
            logger.info(f"L3T_L4T_JET exit code: {exit_code}")
        except ECOSTRESSExitCodeException as e:
            logger.exception(e)
            logger.error(f"ECOSTRESS exit code: {e.exit_code}")
            continue
        except Exception as e:
            logger.exception(e)
            logger.error("exception not handled by ECOSTRESS exit codes")

            if halt_with_unhandled_exceptions:
                break
            else:
                continue

    # scene_df = find_ECOSTRESS_C1_geometry(
    #     target_geometry=geometry,
    #     start_date_UTC=start_date_UTC,
    #     end_date_UTC=end_date_UTC,
    #     start_date_solar=start_date_solar,
    #     end_date_solar=end_date_solar,
    #     hour_min=hour_min,
    #     hour_max=hour_max,
    #     max_cloud_percent=max_cloud_percent,
    #     results_filename=results_filename,
    #     bboxes_filenames=bboxes_filenames
    # )

    # logger.info(
    #     f"processing ECOSTRESS Collection 2 tile {cl.place(tile)} from {cl.time(start_date_UTC)} to {cl.time(end_date_UTC)}")

    # for i, (
    #         time_UTC, time_solar, hour, orbit, scene, cloud_percent, L1B_GEO_filename, L1B_RAD_filename,
    #         L2_LSTE_C1_filename,
    #         L2_CLOUD_C1_filename, area, geometry) in scene_df.iterrows():
    #     logger.info(
    #         f"processing ECOSTRESS Collection 2 tile {cl.place(tile)} from orbit {cl.val(orbit)} scene {cl.val(scene)} at time {cl.time(time_UTC)} UTC {cl.time(time_solar)} solar")
    #     logger.info(f"L1B GEO file: {cl.file(L1B_GEO_filename)}")
    #     logger.info(f"L1B RAD file: {cl.file(L1B_RAD_filename)}")
    #     logger.info(f"L2 LSTE C1 file: {cl.file(L2_LSTE_C1_filename)}")
    #     logger.info(f"L2 CLOUD C1 file: {cl.file(L2_CLOUD_C1_filename)}")

    #     timestamp = f"{time_UTC:%Y%m%dT%H%M%S}"
    #     run_ID = f"ECOv002_{orbit:05d}_{scene:03d}_{timestamp}"

    #     run_directory = join(main_directory, "runs", run_ID)
    #     input_directory = join(run_directory, "input")
    #     makedirs(input_directory, exist_ok=True)
    #     L1B_GEO_input_filename = join(input_directory, basename(L1B_GEO_filename))
    #     copyfile(L1B_GEO_filename, L1B_GEO_input_filename)
    #     L1B_RAD_input_filename = join(input_directory, basename(L1B_RAD_filename))
    #     copyfile(L1B_RAD_filename, L1B_RAD_input_filename)

    #     if output_directory is None:
    #         output_directory = join(main_directory, "output")

    #     makedirs(output_directory, exist_ok=True)

    #     time_UTC = parser.parse(str(time_UTC))
    #     date_UTC = time_UTC.date()

    #     L2_LSTE_output = join(output_directory, "L2_LSTE", date_UTC.strftime("%Y-%m-%d"))
    #     L1_L2_RAD_LSTE_output = join(output_directory, "L1_L2_RAD_LSTE", date_UTC.strftime("%Y-%m-%d"))
    #     L2T_STARS_output = join(output_directory, "L2T_STARS", date_UTC.strftime("%Y-%m-%d"))
    #     L3T_L4T_JET_output = join(output_directory, "L3T_L4T_JET", date_UTC.strftime("%Y-%m-%d"))
    #     L3T_L4T_ALEXI_output = join(output_directory, "L3T_L4T_ALEXI", date_UTC.strftime("%Y-%m-%d"))

    #     # TODO search for previous L2 LSTE and L2 CLOUD

    #     try:
    #         L2_LSTE_granule = L2_LSTE_from_L1B(
    #             L1B_RAD_filename=L1B_RAD_filename,
    #             L1B_GEO_filename=L1B_GEO_filename,
    #             working_directory=run_directory,
    #             native_input_directory=input_directory,
    #             internal_input_directory="/input",
    #             native_output_directory=L2_LSTE_output,
    #             internal_output_directory="/output"
    #         )

    #         L1B_GEO_filename = L2_LSTE_granule.L1B_GEO_filename
    #         L2_LSTE_filename = L2_LSTE_granule.L2_LSTE_filename
    #         L2_CLOUD_filename = L2_LSTE_granule.L2_CLOUD_filename
    #     except Exception as e:
    #         logger.exception(e)
    #         logger.warning(f"L2 LSTE not produced for tile {tile} at time {time_UTC} UTC")
    #         continue

    #     logger.info("preparing to run L1_L2_RAD_LSTE")
    #     logger.info(f"L1B GEO file: {cl.file(L1B_GEO_filename)}")
    #     logger.info(f"L2 LSTE file: {cl.file(L2_LSTE_filename)}")
    #     logger.info(f"L2 CLOUD file: {cl.file(L2_CLOUD_filename)}")

    #     try:
    #         L1_L2_RAD_LSTE_runconfig_filename = generate_L1_L2_RAD_LSTE_runconfig(
    #             L1B_GEO_filename=L1B_GEO_filename,
    #             L1B_RAD_filename=L1B_RAD_filename,
    #             L2_LSTE_filename=L2_LSTE_filename,
    #             L2_CLOUD_filename=L2_CLOUD_filename,
    #             output_directory=L1_L2_RAD_LSTE_output,
    #             working_directory=run_directory,
    #             build=build
    #         )

    #         L1_L2_RAD_LSTE_runconfig = L2GL2TRADLSTEConfig(L1_L2_RAD_LSTE_runconfig_filename)
    #         L2G_LSTE_filename = L1_L2_RAD_LSTE_runconfig.L2G_LSTE_filename
    #         L2G_CLOUD_filename = L1_L2_RAD_LSTE_runconfig.L2G_CLOUD_filename
    #         exit_code = L1_L2_RAD_LSTE(runconfig_filename=L1_L2_RAD_LSTE_runconfig_filename, tiles=[tile])
    #         logger.info(f"L1_L2_RAD_LSTE exit code: {exit_code}")
    #     except ECOSTRESSExitCodeException as e:
    #         logger.exception(e)
    #         logger.error(f"ECOSTRESS exit code: {e.exit_code}")
    #         continue
    #     except Exception as e:
    #         logger.exception(e)
    #         logger.error("exception not handled by ECOSTRESS exit codes")

    #         if halt_with_unhandled_exceptions:
    #             break
    #         else:
    #             continue

    #     pattern = join(abspath(expanduser(L1_L2_RAD_LSTE_output)),
    #                    f"*_L2T_LSTE_{int(orbit):05d}_{int(scene):03d}_{tile}_*.zip")
    #     logger.info(f"searching pattern: {cl.val(pattern)}")
    #     L2T_LSTE_filenames = sorted(glob(pattern))

    #     if len(L2T_LSTE_filenames) == 0:
    #         logger.warning(f"L2T LSTE not produced for tile {tile} at time {time_UTC} UTC")
    #         continue

    #     L2T_LSTE_filename = L2T_LSTE_filenames[0]

    #     logger.info(f"L2T LSTE: {cl.file(L2T_LSTE_filename)}")
    #     L2T_LSTE_granule = L2TLSTE(L2T_LSTE_filename)

    #     L2T_STARS_indices_directory = join(main_directory, "L2T_STARS_INDICES")
    #     L2T_STARS_model_directory = join(main_directory, "L2T_STARS_MODEL")
    #     L2T_STARS_filename = None

    #     try:
    #         L2T_STARS_runconfig_filename = generate_L2T_STARS_runconfig(
    #             orbit=orbit,
    #             scene=scene,
    #             tile=tile,
    #             L2T_LSTE_filename=L2T_LSTE_filename,
    #             working_directory=run_directory,
    #             sources_directory=L2T_STARS_sources_directory,
    #             indices_directory=L2T_STARS_indices_directory,
    #             model_directory=L2T_STARS_model_directory,
    #             output_directory=L2T_STARS_output
    #         )

    #         L2T_STARS_config = L2TSTARSConfig(L2T_STARS_runconfig_filename)
    #         L2T_STARS_filename = L2T_STARS_config.L2T_STARS_zip_filename
    #         exit_code = L2T_STARS(runconfig_filename=L2T_STARS_runconfig_filename, use_VNP43NRT=False)
    #         logger.info(f"L2T_STARS exit code: {exit_code}")
    #     except ECOSTRESSExitCodeException as e:
    #         logger.exception(e)
    #         logger.error(f"ECOSTRESS exit code: {e.exit_code}")
    #         continue
    #     except Exception as e:
    #         logger.exception(e)
    #         logger.error("exception not handled by ECOSTRESS exit codes")

    #         if halt_with_unhandled_exceptions:
    #             break
    #         else:
    #             continue

    #     if not exists(L2T_LSTE_filename):
    #         logger.warning(f"L2T LSTE file not found: {L2T_LSTE_filename}")
    #         continue

    #     if L2T_STARS_filename is None or not exists(L2T_STARS_filename):
    #         logger.warning(f"L2T STARS file not found: {L2T_STARS_filename}")
    #         continue

    #     L1B_GEO_granule = L1BGEO(L1B_GEO_filename)
    #     average_SZA_degrees = float(L1B_GEO_granule.product_metadata["AverageSolarZenith"])

    #     if average_SZA_degrees > SZA_cutoff:
    #         logger.warning(
    #             f"skipping L3/L4 with SZA {average_SZA_degrees} degrees for orbit {orbit} scene {scene} tile {tile} at {time_UTC} UTC {time_solar} solar")
    #         continue

    #     try:
    #         L3T_L4T_JET_runconfig_filename = generate_L3T_L4T_JET_runconfig(
    #             orbit=orbit,
    #             scene=scene,
    #             tile=tile,
    #             L2T_LSTE_filename=L2T_LSTE_filename,
    #             L2T_STARS_filename=L2T_STARS_filename,
    #             working_directory=run_directory,
    #             sources_directory=L3T_L4T_JET_sources_directory,
    #             static_directory=static_directory,
    #             output_directory=L3T_L4T_JET_output,
    #             SRTM_directory=SRTM_directory
    #         )

    #         L3T_L4T_JET_config = L3TL4TJETConfig(L3T_L4T_JET_runconfig_filename)
    #         L3T_L4T_JET_sources_directory = L3T_L4T_JET_config.output_directory
    #         logger.info(f"L3T L4T directory: {L3T_L4T_JET_sources_directory}")
    #         exit_code = L3T_L4T_JET(runconfig_filename=L3T_L4T_JET_runconfig_filename)
    #         logger.info(f"L3T_L4T_JET exit code: {exit_code}")
    #     except ECOSTRESSExitCodeException as e:
    #         logger.exception(e)
    #         logger.error(f"ECOSTRESS exit code: {e.exit_code}")
    #         continue
    #     except Exception as e:
    #         logger.exception(e)
    #         logger.error("exception not handled by ECOSTRESS exit codes")

    #         if halt_with_unhandled_exceptions:
    #             break
    #         else:
    #             continue

    #     try:
    #         L3T_L4T_ALEXI_runconfig_filename = generate_L3T_L4T_ALEXI_runconfig(
    #             orbit=orbit,
    #             scene=scene,
    #             tile=tile,
    #             L2T_LSTE_filename=L2T_LSTE_filename,
    #             L2T_STARS_filename=L2T_STARS_filename,
    #             working_directory=run_directory,
    #             sources_directory=ALEXI_directory,
    #             output_directory=L3T_L4T_ALEXI_output,
    #             static_directory=static_directory
    #         )

    #         L3T_L4T_ALEXI_config = L3TL4TALEXIConfig(L3T_L4T_ALEXI_runconfig_filename)
    #         L3T_L4T_ALEXI_directory = L3T_L4T_ALEXI_config.output_directory
    #         logger.info(f"L3T L4T disALEXI directory: {L3T_L4T_ALEXI_directory}")
    #         exit_code = L3T_L4T_ALEXI(runconfig_filename=L3T_L4T_ALEXI_runconfig_filename)
    #         logger.info(f"L3T_L4T_ALEXI exit code: {exit_code}")
    #     except ECOSTRESSExitCodeException as e:
    #         logger.exception(e)
    #         logger.error(f"ECOSTRESS exit code: {e.exit_code}")
    #         continue
    #     except Exception as e:
    #         logger.exception(e)
    #         logger.error("exception not handled by ECOSTRESS exit codes")

    #         if halt_with_unhandled_exceptions:
    #             break
    #         else:
    #             continue


def main(argv=sys.argv):
    tile = argv[1]

    if "--start-date-UTC" in argv:
        start_date_UTC = parser.parse(argv[argv.index("--start-date-UTC") + 1]).date()
    elif "--start" in argv:
        start_date_UTC = parser.parse(argv[argv.index("--start") + 1]).date()
    else:
        start_date_UTC = None

    if "--end-date-UTC" in argv:
        end_date_UTC = parser.parse(argv[argv.index("--end-date-UTC") + 1]).date()
    elif "--end" in argv:
        end_date_UTC = parser.parse(argv[argv.index("--end") + 1]).date()
    else:
        end_date_UTC = None

    if "--date" in argv:
        date_UTC = parser.parse(argv[argv.index("--date") + 1]).date()
        start_date_UTC = date_UTC
        end_date_UTC = date_UTC

    if "--max-cloud" in argv:
        max_cloud_percent = float(argv[argv.index("--max-cloud") + 1])
    else:
        max_cloud_percent = None

    if "--results" in argv:
        results_filename = str(argv[argv.index("--results") + 1])
    else:
        results_filename = None

    if "--main" in argv:
        main_directory = str(argv[argv.index("--main") + 1])
    else:
        main_directory = None

    if "--static" in argv:
        static_directory = str(argv[argv.index("--static") + 1])
    else:
        static_directory = None

    if "--SRTM" in argv:
        SRTM_directory = str(argv[argv.index("--SRTM") + 1])
    else:
        SRTM_directory = None

    if "--JET" in argv:
        L2T_STARS_sources_directory = str(argv[argv.index("--JET") + 1])
    else:
        L2T_STARS_sources_directory = None

    if "--STARS" in argv:
        L3T_L4T_JET_sources_directory = str(argv[argv.index("--STARS") + 1])
    else:
        L3T_L4T_JET_sources_directory = None

    if "--ALEXI" in argv:
        ALEXI_directory = str(argv[argv.index("--ALEXI") + 1])
    else:
        ALEXI_directory = None

    if "--output" in argv:
        output_directory = str(argv[argv.index("--output") + 1])
    else:
        output_directory = None

    reprocess_JET(
        tile=tile,
        start_date_UTC=start_date_UTC,
        end_date_UTC=end_date_UTC,
        max_cloud_percent=max_cloud_percent,
        results_filename=results_filename,
        main_directory=main_directory,
        static_directory=static_directory,
        SRTM_directory=SRTM_directory,
        L2T_STARS_sources_directory=L2T_STARS_sources_directory,
        L3T_L4T_JET_sources_directory=L3T_L4T_JET_sources_directory,
        ALEXI_directory=ALEXI_directory,
        output_directory=output_directory
    )


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
