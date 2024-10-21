import logging

import numpy as np

from ECOSTRESS.L3_ET_PTJPL import L3ETPTJPL, L3GETPTJPL


def L3T_ET_PTJPL_simulation(
        L1B_GEO_filename: str,
        L2_CLOUD_filename: str,
        L3_ET_PTJPL_filename: str,
        gridded_products_directory: str = "gridded",
        tiled_products_directory: str = "tiled"):
    logger = logging.getLogger(__name__)

    L3_PTJPL_granule = L3ETPTJPL(
        L3_ET_PTJPL_filename=L3_ET_PTJPL_filename,
        L2_CLOUD_filename=L2_CLOUD_filename,
        L1B_GEO_filename=L1B_GEO_filename
    )

    L3G_PTJPL_granule = L3GETPTJPL.from_swath(
        L3_PTJPL_granule,
        output_directory=gridded_products_directory
    )

    L3G_PTJPL_storage = L3G_PTJPL_granule.storage
    L3T_PTJPL_tiles = L3G_PTJPL_granule.to_tiles(output_directory=tiled_products_directory)
    L3T_PTJPL_tile_sizes = [tile.storage for tile in L3T_PTJPL_tiles]
    L3T_PTJPL_scene_storage = np.nansum(L3T_PTJPL_tile_sizes)
    L3T_PTJPL_tile_storage = np.nanmax(L3T_PTJPL_tile_sizes)

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
#         L4_PTJPL_ESI_filename=files["L4_ESI_PT-JPL"],
#         L4_disALEXI_ESI_filename=files["L4_ESI_ALEXI"],
#         L4_PTJPL_WUE_filename=files["L4_WUE"]
#     )
