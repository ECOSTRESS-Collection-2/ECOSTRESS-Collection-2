import logging
from datetime import datetime, date
from os.path import exists
from typing import List
import h5py
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon
from dateutil import parser
import colored_logging as cl
import rasters as rt
from BRDF import bidirectional_reflectance
from BRDF.SZA import calculate_SZA
from GEOS5FP.GEOS5FP import GEOS5FP
from MODLAND import find_MODLAND_tiles
from MODLAND.indices import generate_MODLAND_grid, parsehv
from rasters import Raster, RasterGrid, RasterGeometry
from .VIIRSDownloader import VIIRSDownloaderAlbedo
from .VIIRSDataPool import VIIRSDataPool, VIIRSGranule

NDVI_COLORMAP = "jet_r"
ALBEDO_COLORMAP = "gray"

logger = logging.getLogger(__name__)


class VIIRSUnavailableError(Exception):
    pass


class VNP43MA3Granule(VIIRSGranule):
    def __init__(
            self,
            filename: str,
            working_directory: str = None,
            products_directory: str = None,
            GEOS5FP_connection: GEOS5FP = None,
            GEOS5FP_download: str = None,
            GEOS5FP_products: str = None):
        super(VNP43MA3Granule, self).__init__(
            filename=filename,
            working_directory=working_directory,
            products_directory=products_directory
        )

        if GEOS5FP_connection is None:
            GEOS5FP_connection = GEOS5FP(
                working_directory=working_directory,
                download_directory=GEOS5FP_download,
                products_directory=GEOS5FP_products
            )

        self.GEOS5FP = GEOS5FP_connection

    def BSA(
            self,
            band: int,
            geometry: RasterGeometry = None,
            save_data: bool = True,
            include_preview: bool = True,
            apply_QA: bool = True,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"BSA_M{band}")

        if product_filename is not None and exists(product_filename):
            logger.info(f"loading VNP43MA3 BSA M{band}: {cl.file(product_filename)}")
            image = Raster.open(product_filename)
        else:
            image = self.dataset(
                filename=self.filename,
                dataset_name=f"HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/Albedo_BSA_M{int(band)}",
                fill_value=32767,
                scale_factor=0.001
            )

            if apply_QA:
                QA = self.QA(
                    band=band,
                    geometry=geometry,
                    save_data=save_data,
                    include_preview=include_preview
                )

                image = rt.where(QA == 0, image, np.nan)

        if save_data and not exists(product_filename):
            logger.info(f"writing VNP43MA3 BSA M{band}: {cl.file(product_filename)}")
            image.to_geotiff(product_filename, include_preview=include_preview)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    def WSA(
            self,
            band: int,
            geometry: RasterGeometry = None,
            save_data: bool = True,
            include_preview: bool = True,
            apply_QA: bool = True,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"WSA_M{band}")

        if product_filename is not None and exists(product_filename):
            logger.info(f"loading VNP43MA3 WSA M{band}: {cl.file(product_filename)}")
            image = Raster.open(product_filename)
        else:
            image = self.dataset(
                filename=self.filename,
                dataset_name=f"HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/Albedo_WSA_M{int(band)}",
                fill_value=32767,
                scale_factor=0.001
            )

            if apply_QA:
                QA = self.QA(
                    band=band,
                    geometry=geometry,
                    save_data=save_data,
                    include_preview=include_preview
                )

                image = rt.where(QA == 0, image, np.nan)

        if save_data and not exists(product_filename):
            logger.info(f"writing VNP43MA3 WSA M{band}: {cl.file(product_filename)}")
            image.to_geotiff(product_filename, include_preview=include_preview)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    def QA(
            self,
            band: int,
            geometry: RasterGeometry = None,
            save_data: bool = True,
            include_preview: bool = True,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"VNP43MA3_QA_M{band}")

        if product_filename is not None and exists(product_filename):
            logger.info(f"loading VNP43MA3 QA M{band}: {cl.file(product_filename)}")
            image = Raster.open(product_filename)
        else:
            dataset_name = f"HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M{int(band)}"

            with h5py.File(self.filename, "r") as f:
                image = np.array(f[dataset_name])
                h, v = self.hv
                grid = generate_MODLAND_grid(h, v, image.shape[0])
                logger.info(f"opening file: {cl.file(self.filename)}")
                logger.info(f"loading {cl.val(dataset_name)} at {cl.val(f'{grid.cell_size:0.2f} m')} resolution")
                image = Raster(image, geometry=grid)

        if save_data and not exists(product_filename):
            logger.info(f"writing VNP43MA3 QA M{band}: {cl.file(product_filename)}")
            image.to_geotiff(product_filename, include_preview=include_preview)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    @property
    def geometry(self) -> RasterGrid:
        return generate_MODLAND_grid(*parsehv(self.tile), 1200)

    def get_albedo(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = True,
            include_preview: bool = True,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename("albedo")

        if product_filename is not None and exists(product_filename):
            logger.info(f"loading VNP43MA3 albedo: {cl.file(product_filename)}")
            image = Raster.open(product_filename)
        else:
            date_UTC = self.date_UTC
            doy = date_UTC.timetuple().tm_yday
            SZA = calculate_SZA(doy, 10.5, self.geometry)
            time_UTC = datetime(date_UTC.year, date_UTC.month, date_UTC.day, 10, 30)
            AOT = self.GEOS5FP.AOT(time_UTC=time_UTC, geometry=self.geometry, resampling="cubic")

            b = {}

            for m in (1, 2, 3, 4, 5, 7, 8, 10, 11):
                WSA = self.WSA(m)
                BSA = self.BSA(m)

                band_albedo = bidirectional_reflectance(
                    white_sky_albedo=WSA,
                    black_sky_albedo=BSA,
                    SZA=SZA,
                    AOT=AOT
                )

                b[m] = band_albedo

            image = 0.2418 * b[1] \
                    - 0.201 * b[2] \
                    + 0.2093 * b[3] \
                    + 0.1146 * b[4] \
                    + 0.1348 * b[5] \
                    + 0.2251 * b[7] \
                    + 0.1123 * b[8] \
                    + 0.0860 * b[10] \
                    + 0.0803 * b[11] \
                    - 0.0131

        if save_data and not exists(product_filename):
            logger.info(f"writing VNP43MA3 albedo: {cl.file(product_filename)}")
            image.to_geotiff(product_filename, include_preview=include_preview)

        if geometry is not None:
            image = image.to_geometry(geometry)

        image.cmap = ALBEDO_COLORMAP

        return image

    albedo = property(get_albedo)


class VNP43MA3(VIIRSDataPool, VIIRSDownloaderAlbedo):
    DEFAULT_DOWNLOAD_DIRECTORY = "VNP43MA3_download"
    DEFAULT_PRODUCTS_DIRECTORY = "VNP43MA3_products"
    DEFAULT_MOSAIC_DIRECTORY = "VNP43MA3_mosaics"

    def search(
            self,
            start_date: date or datetime or str,
            end_date: date or datetime or str = None,
            build: str = None,
            tiles: List[str] or str = None,
            target_geometry: Point or Polygon or RasterGrid = None,
            *args,
            **kwargs) -> pd.DataFrame:
        return super(VNP43MA3, self).search(
            product="VNP43MA3",
            start_date=start_date,
            end_date=end_date,
            build=build,
            tiles=tiles,
            target_geometry=target_geometry,
            *args,
            **kwargs
        )

    def granule(
            self,
            date_UTC: date or str,
            tile: str,
            download_location: str = None,
            build: str = None) -> VNP43MA3Granule:
        listing = self.search(
            start_date=date_UTC,
            end_date=date_UTC,
            build=build,
            tiles=[tile]
        )

        if len(listing) > 0:
            URL = listing.iloc[0].URL

        filename = super(VNP43MA3, self).download_URL(
            URL=URL,
            download_location=download_location
        )

        granule = VNP43MA3Granule(
            filename=filename,
            products_directory=self.products_directory
        )

        return granule

    def albedo(
            self,
            date_UTC: date or str,
            geometry: RasterGeometry,
            filename: str = None,
            resampling: str = None) -> Raster:
        if isinstance(date_UTC, str):
            date_UTC = parser.parse(date_UTC).date()

        if filename is not None and exists(filename):
            return Raster.open(filename, cmap=ALBEDO_COLORMAP)

        # if resampling is None:
        #     resampling = self.resampling

        tiles = sorted(find_MODLAND_tiles(geometry.boundary_latlon.geometry))
        albedo = None

        for tile in tiles:
            granule = self.granule(date_UTC=date_UTC, tile=tile)
            granule_albedo = granule.albedo
            source_cell_size = granule_albedo.geometry.cell_size
            dest_cell_size = geometry.cell_size
            logger.info(f"projecting VIIRS albedo from {cl.val(f'{source_cell_size} m')} to {cl.val(f'{dest_cell_size} m')}")
            projected_albedo = granule_albedo.to_geometry(geometry, resampling=resampling)

            if albedo is None:
                albedo = projected_albedo
            else:
                albedo = rt.where(np.isnan(albedo), projected_albedo, albedo)

        albedo.cmap = ALBEDO_COLORMAP

        if filename is not None:
            logger.info(f"writing albedo mosaic: {cl.file(filename)}")
            albedo.to_geotiff(filename)

        return albedo
