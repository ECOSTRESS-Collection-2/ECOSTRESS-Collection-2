import logging

import numpy as np

from ECOSTRESS.L4_WUE import L4WUE, L4GWUE


def L4T_WUE_simulation(
        L1B_GEO_filename: str,
        L2_CLOUD_filename: str,
        L4_WUE_filename: str,
        gridded_products_directory: str = "gridded",
        tiled_products_directory: str = "tiled"):
    logger = logging.getLogger(__name__)

    L4_WUE_granule = L4WUE(
        L4_WUE_filename=L4_WUE_filename,
        L2_CLOUD_filename=L2_CLOUD_filename,
        L1B_GEO_filename=L1B_GEO_filename
    )

    L4G_WUE_granule = L4GWUE.from_swath(
        L4_WUE_granule,
        output_directory=gridded_products_directory
    )

    L4G_WUE_storage = L4G_WUE_granule.storage
    L4T_WUE_tiles = L4G_WUE_granule.to_tiles(output_directory=tiled_products_directory)
    L4T_WUE_tile_sizes = [tile.storage for tile in L4T_WUE_tiles]
    L4T_WUE_scene_storage = np.nansum(L4T_WUE_tile_sizes)
    L4T_WUE_tile_storage = np.nanmax(L4T_WUE_tile_sizes)

# def main(argv=sys.argv):
#     configure_logger()
#     orbit = int(argv[1])
#     scene = int(argv[2])
#     files = find_ECOSTRESS(orbit, scene)["0601"]
#
#     gridded_tiled_simulation(
#         L1B_GEO_filename=files["L1B_GEO"],
#         L2_CLOUD_filename=files["L2_CLOUD"],
#         L2_LSTE_filename=files["L2_LSTE"],
#         L3_PTJPL_filename=files["L3_ET_PT-JPL"],
#         L3_disALEXI_filename=files["L3_ET_ALEXI"],
#         L4_PTJPL_WUE_filename=files["L4_WUE_PT-JPL"],
#         L4_disALEXI_WUE_filename=files["L4_WUE_ALEXI"],
#         L4_PTJPL_WUE_filename=files["L4_WUE"]
#     )
