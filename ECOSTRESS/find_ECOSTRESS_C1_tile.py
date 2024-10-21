import sys

from dateutil import parser

from ECOSTRESS.find_ECOSTRESS_C1_geometry import find_ECOSTRESS_C1_geometry
from sentinel_tile_grid import sentinel_tile_grid


def main(argv=sys.argv):
    tile = argv[1]

    geometry = sentinel_tile_grid.grid(tile)

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

    if "--max-cloud" in argv:
        max_cloud_percent = float(argv[argv.index("--max-cloud") + 1])
    else:
        max_cloud_percent = None

    if "--results" in argv:
        results_filename = parser.parse(argv[argv.index("--results") + 1])
    else:
        results_filename = f"{tile}_{start_date_UTC}~{end_date_UTC}.geojson"

    find_ECOSTRESS_C1_geometry(
        target_geometry=geometry,
        start_date_UTC=start_date_UTC,
        end_date_UTC=end_date_UTC,
        max_cloud_percent=max_cloud_percent,
        results_filename=results_filename
    )


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
