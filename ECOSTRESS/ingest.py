from ECOSTRESS.ECOSTRESS_swath_granule import ECOSTRESSSwathGranule
from rasters import Raster


def read_ECOSTRESS_dataset(
        product_filename: str,
        dataset_name: str,
        geolocation_filename: str,
        cloud_filename: str = None) -> Raster:
    granule = ECOSTRESSSwathGranule(
        product_filename=product_filename,
        L2_CLOUD_filename=cloud_filename,
        L1B_GEO_filename=geolocation_filename
    )

    data = granule.data(dataset_name)

    return data
