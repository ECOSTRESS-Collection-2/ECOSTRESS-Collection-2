import logging
import warnings
from datetime import datetime, date
from os import remove
from os.path import exists, join, abspath, expanduser
import re
from pathlib import Path
import tempfile
from typing import List, Union

import earthaccess
import h5py
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from daterange import get_date
from dateutil import parser
from skimage.transform import resize

import colored_logging as cl
import ecostress_cmr
from ECOSTRESS.exit_codes import DownloadFailed
import rasters
from LPDAAC.LPDAACDataPool import RETRIES
from MODLAND.indices import generate_MODLAND_grid
from rasters import Raster, RasterGrid, RasterGeometry
from rasters import Point, Polygon
from .VIIRSDataPool import VIIRSGranule

NDVI_COLORMAP = LinearSegmentedColormap.from_list(
    name="NDVI",
    colors=[
        "#0000ff",
        "#000000",
        "#745d1a",
        "#e1dea2",
        "#45ff01",
        "#325e32"
    ]
)

ALBEDO_COLORMAP = "gray"

logger = logging.getLogger(__name__)

class VIIRSUnavailableError(Exception):
    pass


class VNP09GAGranule(VIIRSGranule):
    CLOUD_DATASET_NAME = "HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SurfReflect_QF1_1"

    def get_cloud_mask(self, target_shape: tuple = None) -> Raster:
        h, v = self.hv

        if self._cloud_mask is None:
            with h5py.File(self.filename, "r") as f:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")
                    QF1 = np.array(f[self.CLOUD_DATASET_NAME])

                cloud_levels = (QF1 >> 2) & 3
                cloud_mask = cloud_levels > 0
                self._cloud_mask = cloud_mask
        else:
            cloud_mask = self._cloud_mask

        if target_shape is not None:
            cloud_mask = resize(cloud_mask, target_shape, order=0).astype(bool)
            shape = target_shape
        else:
            shape = cloud_mask.shape

        geometry = generate_MODLAND_grid(h, v, shape[0])
        cloud_mask = Raster(cloud_mask, geometry=geometry)

        return cloud_mask

    cloud_mask = property(get_cloud_mask)

    def dataset(
            self,
            filename: str,
            dataset_name: str,
            scale_factor: float,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            resampling: str = None) -> Raster:

        with h5py.File(filename, "r") as f:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                dataset = f[dataset_name]
                DN = np.array(dataset)

            if "_FillValue" in dataset.attrs:
                fill_value = dataset.attrs["_FillValue"]
            else:
                fill_value = dataset.attrs["_Fillvalue"]

            h, v = self.hv
            grid = generate_MODLAND_grid(h, v, DN.shape[0])
            logger.info(f"opening VIIRS file: {cl.file(self.filename)}")
            logger.info(f"loading {cl.val(dataset_name)} at {cl.val(f'{grid.cell_size:0.2f} m')} resolution")
            DN = np.where(DN == fill_value, np.nan, DN)
            DN = Raster(DN, geometry=grid)

        data = DN * scale_factor

        if apply_cloud_mask:
            if cloud_mask is None:
                cloud_mask = self.get_cloud_mask(target_shape=DN.shape)

            data = rasters.where(cloud_mask, np.nan, data)

        if geometry is not None:
            data = data.to_geometry(geometry, resampling=resampling)

        return data

    @property
    def geometry_M(self) -> RasterGrid:
        return generate_MODLAND_grid(*self.hv, 1200)

    @property
    def geometry_I(self) -> RasterGrid:
        return generate_MODLAND_grid(*self.hv, 2400)

    def geometry(self, band: str) -> RasterGrid:
        try:
            band_letter = band[0]
        except Exception as e:
            raise ValueError(f"invalid band: {band}")

        if band_letter == "I":
            return self.geometry_I
        elif band_letter == "M":
            return self.geometry_M
        else:
            raise ValueError(f"invalid band: {band}")

    def get_sensor_zenith_M(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"sensor_zenith_M")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS sensor zenith: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SensorZenith_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank sensor zenith image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS M-band sensor zenith: {cl.file(product_filename)} {cl.val(image.shape)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    sensor_zenith_M = property(get_sensor_zenith_M)

    def get_sensor_zenith_I(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"sensor_zenith_I")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS I-band sensor zenith: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            h, v = self.hv
            grid_I = generate_MODLAND_grid(h, v, 2400)

            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SensorZenith_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False,
                geometry=grid_I,
                resampling="cubic"
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank sensor zenith image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS sensor zenith: {cl.file(product_filename)} {cl.val(image.shape)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    sensor_zenith_I = property(get_sensor_zenith_I)

    def sensor_zenith(
            self,
            band: str,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        try:
            band_letter = band[0]
        except Exception as e:
            raise ValueError(f"invalid band: {band}")

        if band_letter == "I":
            return self.get_sensor_zenith_I(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        elif band_letter == "M":
            return self.get_sensor_zenith_M(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        else:
            raise ValueError(f"invalid band: {band}")

    def get_sensor_azimuth_M(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"sensor_azimuth_M")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS M-band sensor azimuth: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SensorAzimuth_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank sensor azimuth image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS sensor azimuth: {cl.file(product_filename)} {cl.val(image.shape)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    sensor_azimuth_M = property(get_sensor_azimuth_M)

    def get_sensor_azimuth_I(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"sensor_azimuth_I")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS I-band sensor azimuth: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            h, v = self.hv
            grid_I = generate_MODLAND_grid(h, v, 2400)

            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SensorAzimuth_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False,
                geometry=grid_I,
                resampling="cubic"
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank sensor azimuth image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS sensor azimuth: {cl.file(product_filename)} {cl.val(image.shape)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    sensor_azimuth_I = property(get_sensor_azimuth_I)

    def sensor_azimuth(
            self,
            band: str,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        try:
            band_letter = band[0]
        except Exception as e:
            raise ValueError(f"invalid band: {band}")

        if band_letter == "I":
            return self.get_sensor_azimuth_I(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        elif band_letter == "M":
            return self.get_sensor_azimuth_M(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        else:
            raise ValueError(f"invalid band: {band}")

    def get_solar_zenith_M(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"solar_zenith_M")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS M-band solar zenith: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SolarZenith_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank solar zenith image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS solar zenith: {cl.file(product_filename)} {cl.val(image.shape)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    solar_zenith_M = property(get_solar_zenith_M)

    def get_solar_zenith_I(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename("solar_zenith_I")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS I-band solar zenith: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            h, v = self.hv
            grid_I = generate_MODLAND_grid(h, v, 2400)

            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SolarZenith_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False,
                geometry=grid_I,
                resampling="cubic"
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank solar zenith image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS solar zenith: {cl.file(product_filename)} {cl.val(image.shape)}")

            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    solar_zenith_I = property(get_solar_zenith_I)

    def solar_zenith(
            self,
            band: str,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        try:
            band_letter = band[0]
        except Exception as e:
            raise ValueError(f"invalid band: {band}")

        if band_letter == "I":
            return self.get_solar_zenith_I(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        elif band_letter == "M":
            return self.get_solar_zenith_M(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        else:
            raise ValueError(f"invalid band: {band}")

    def get_solar_azimuth_M(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename("solar_azimuth_M")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS M-band solar azimuth: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SolarAzimuth_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank solar azimuth image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS solar azimuth: {cl.file(product_filename)} {cl.val(image.shape)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    solar_azimuth_M = property(get_solar_azimuth_M)

    def get_solar_azimuth_I(
            self,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename("solar_azimuth_I")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS I-band solar azimuth: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            h, v = self.hv
            grid_I = generate_MODLAND_grid(h, v, 2400)

            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SolarAzimuth_1",
                0.01,
                cloud_mask=None,
                apply_cloud_mask=False,
                geometry=grid_I,
                resampling="cubic"
            )

        if np.all(np.isnan(image)):
            raise ValueError("blank solar azimuth image")

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS solar azimuth: {cl.file(product_filename)} {cl.val(image.shape)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    solar_azimuth_I = property(get_solar_azimuth_I)

    def solar_azimuth(
            self,
            band: str,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        try:
            band_letter = band[0]
        except Exception as e:
            raise ValueError(f"invalid band: {band}")

        if band_letter == "I":
            return self.get_solar_azimuth_I(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        elif band_letter == "M":
            return self.get_solar_azimuth_M(
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        else:
            raise ValueError(f"invalid band: {band}")

    def get_M_band(
            self,
            band: int,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"M{band}")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS M-band {band} surface reflectance: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_1km_2D/Data Fields/SurfReflect_M{int(band)}_1",
                0.0001,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask
            )

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS M{band}: {cl.file(product_filename)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    def get_I_band(
            self,
            band: int,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename(f"I{band}")

        image = None

        if product_filename is not None and exists(product_filename):
            try:
                logger.info(f"loading VIIRS I-band {band} surface reflectance: {cl.file(product_filename)}")
                image = Raster.open(product_filename)
            except Exception as e:
                logger.exception(e)
                logger.warning(f"removing corrupted file: {product_filename}")
                remove(product_filename)
                image = None

        if image is None:
            image = self.dataset(
                self.filename,
                f"HDFEOS/GRIDS/VIIRS_Grid_500m_2D/Data Fields/SurfReflect_I{int(band)}_1",
                0.0001,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask
            )

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS I{band}: {cl.file(product_filename)}")
            image.to_geotiff(product_filename)

            if save_preview:
                image.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"), quality=20, remove_XML=True)

        if geometry is not None:
            image = image.to_geometry(geometry)

        return image

    def band(
            self,
            band: str,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        try:
            band_letter = band[0]
            band_number = int(band[1:])
        except Exception as e:
            raise ValueError(f"invalid band: {band}")

        if band_letter == "I":
            return self.get_I_band(
                band=band_number,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        elif band_letter == "M":
            return self.get_M_band(
                band=band_number,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview,
                product_filename=product_filename
            )
        else:
            raise ValueError(f"invalid band: {band}")

    def get_red(
            self,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        return self.get_I_band(
            band=1,
            cloud_mask=cloud_mask,
            apply_cloud_mask=apply_cloud_mask,
            geometry=geometry,
            save_data=save_data,
            save_preview=save_preview,
            product_filename=product_filename
        )

    red = property(get_red)

    def get_NIR(
            self,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        return self.get_I_band(
            band=2,
            cloud_mask=cloud_mask,
            apply_cloud_mask=apply_cloud_mask,
            geometry=geometry,
            save_data=save_data,
            save_preview=save_preview,
            product_filename=product_filename
        )

    NIR = property(get_NIR)

    def get_NDVI(
            self,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename("NDVI")

        if product_filename is not None and exists(product_filename):
            logger.info(f"loading VIIRS NDVI: {cl.file(product_filename)}")
            NDVI = Raster.open(product_filename)
        else:
            red = self.get_red(
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            NIR = self.get_NIR(
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            NDVI = np.clip((NIR - red) / (NIR + red), -1, 1)

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS NDVI: {cl.file(product_filename)}")
            NDVI.to_geotiff(product_filename)

            if save_preview:
                NDVI.percentilecut.to_geojpeg(product_filename.replace(".tif", ".jpeg"))

        if geometry is not None:
            NDVI = NDVI.to_geometry(geometry)

        NDVI.cmap = NDVI_COLORMAP

        return NDVI

    NDVI = property(get_NDVI)

    def get_albedo(
            self,
            cloud_mask: Raster = None,
            apply_cloud_mask: bool = True,
            geometry: RasterGeometry = None,
            save_data: bool = False,
            save_preview: bool = False,
            product_filename: str = None) -> Raster:
        if product_filename is None:
            product_filename = self.product_filename("albedo")

        if product_filename is not None and exists(product_filename):
            logger.info(f"loading VIIRS albedo: {cl.file(product_filename)}")
            albedo = Raster.open(product_filename)
        else:
            b1 = self.get_M_band(
                1,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b2 = self.get_M_band(
                2,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b3 = self.get_M_band(
                3,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b4 = self.get_M_band(
                4,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b5 = self.get_M_band(
                5,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b7 = self.get_M_band(
                7,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b8 = self.get_M_band(
                8,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b10 = self.get_M_band(
                10,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            b11 = self.get_M_band(
                11,
                cloud_mask=cloud_mask,
                apply_cloud_mask=apply_cloud_mask,
                geometry=geometry,
                save_data=save_data,
                save_preview=save_preview
            )

            # https://lpdaac.usgs.gov/documents/194/VNP43_ATBD_V1.pdf
            albedo = 0.2418 * b1 \
                     - 0.201 * b2 \
                     + 0.2093 * b3 \
                     + 0.1146 * b4 \
                     + 0.1348 * b5 \
                     + 0.2251 * b7 \
                     + 0.1123 * b8 \
                     + 0.0860 * b10 \
                     + 0.0803 * b11 \
                     - 0.0131

            albedo = np.clip(albedo, 0, 1)

        if save_data and not exists(product_filename):
            logger.info(f"writing VIIRS albedo: {cl.file(product_filename)}")
            albedo.to_geotiff(product_filename)

        if geometry is not None:
            logger.info(f"projecting VIIRS albedo from {cl.val(albedo.geometry.cell_size)} to {cl.val(geometry.cell_size)}")
            albedo = albedo.to_geometry(geometry)

        albedo.cmap = ALBEDO_COLORMAP

        return albedo

    albedo = property(get_albedo)


VIIRS_CONCEPT = "C2631841556-LPCLOUD"

def earliest_datetime(date_in: Union[date, str]) -> datetime:
    if isinstance(date_in, str):
        datetime_in = parser.parse(date_in)
    else:
        datetime_in = date_in

    date_string = datetime_in.strftime("%Y-%m-%d")
    return parser.parse(f"{date_string}T00:00:00Z")


def latest_datetime(date_in: Union[date, str]) -> datetime:
    if isinstance(date_in, str):
        datetime_in = parser.parse(date_in)
    else:
        datetime_in = date_in

    date_string = datetime_in.strftime("%Y-%m-%d")
    return parser.parse(f"{date_string}T23:59:59Z")


VIIRS_FILENAME_REGEX = re.compile("^VNP09GA\.[^.]+\.([^.]+)\.002\.\d+\.h5$")
def modland_tile_from_filename(filename: str) -> str:
    match = VIIRS_FILENAME_REGEX.match(filename)
    if match is None:
        raise RuntimeError(f"Invalid filename found through VIIRS CMR search: {filename}")

    return match.group(1)


# TODO: Deduplicate between VIIRS and HLS
def VIIRS_CMR_query(
        start_date: Union[date, str],
        end_date: Union[date, str],
        target_geometry: Point or Polygon or RasterGeometry = None,
        tile: str = None,
) -> List[earthaccess.search.DataGranule]:
    """function to search for VIIRS at tile in date range"""
    query = earthaccess.granule_query() \
        .concept_id(VIIRS_CONCEPT) \
        .temporal(earliest_datetime(start_date), latest_datetime(end_date))

    if isinstance(target_geometry, Point):
        query = query.point(target_geometry.x, target_geometry.y)
    if isinstance(target_geometry, Polygon):
        ring = target_geometry.exterior
        if not ring.is_ccw:
            ring = ring.reverse()
        coordinates = ring.coords
        query = query.polygon(coordinates)
    if isinstance(target_geometry, RasterGeometry):
        ring = target_geometry.corner_polygon_latlon.exterior
        if not ring.is_ccw:
            ring = ring.reverse()
        coordinates = ring.coords
        query = query.polygon(coordinates)
    if tile is not None:
        query = query.readable_granule_name(f"*.{tile}.*")

    granules: List[earthaccess.search.DataGranule]
    try:
        granules = query.get()
    except Exception as e:
        raise ecostress_cmr.CMRServerUnreachable(e)
    granules = sorted(granules, key=lambda granule: granule["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"])

    logger.info("Found the following granules for VIIRS 2 using the CMR search:")
    for granule in granules:
        logger.info("  " + cl.file(granule["meta"]["native-id"]))
    logger.info(f"Number of VIIRS 2 granules found using CMR search: {len(granules)}")

    return granules


class VNP09GA:
    DEFAULT_WORKING_DIRECTORY = "."
    DEFAULT_DOWNLOAD_DIRECTORY = "VNP09GA_download"
    DEFAULT_PRODUCTS_DIRECTORY = "VNP09GA_products"
    DEFAULT_MOSAIC_DIRECTORY = "VNP09GA_mosaics"
    DEFAULT_RESAMPLING = "nearest"

    CLOUD_DATASET_NAME = "HDFEOS/GRIDS/VNP_Grid_1km_2D/Data Fields/SurfReflect_QF1_1"

    def __init__(
            self,
            working_directory: str = None,
            download_directory: str = None,
            products_directory: str = None,
            mosaic_directory: str = None,
            resampling: str = None):

        if resampling is None:
            resampling = self.DEFAULT_RESAMPLING

        self.resampling = resampling

        self._granules = pd.DataFrame(["date_UTC", "tile", "granule"])

        if working_directory is None:
            working_directory = self.DEFAULT_WORKING_DIRECTORY

        working_directory = abspath(expanduser(working_directory))

        if download_directory is None:
            download_directory = join(working_directory, self.DEFAULT_DOWNLOAD_DIRECTORY)

        download_directory = abspath(expanduser(download_directory))

        if products_directory is None:
            products_directory = join(working_directory, self.DEFAULT_PRODUCTS_DIRECTORY)

        products_directory = abspath(expanduser(products_directory))

        if mosaic_directory is None:
            mosaic_directory = join(working_directory, self.DEFAULT_MOSAIC_DIRECTORY)

        mosaic_directory = abspath(expanduser(mosaic_directory))

        self.working_directory = working_directory
        self.download_directory = download_directory
        self.products_directory = products_directory
        self.mosaic_directory = mosaic_directory

        self.auth = ecostress_cmr.login()

    def add_granules(self, granules: List[earthaccess.search.DataGranule]):
        data = pd.DataFrame([
            {
                "date_UTC": get_date(granule["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"]),
                "tile": modland_tile_from_filename(Path(granule.data_links()[0]).name),
                "granule": granule,
            }
            for granule in granules
        ])

        self._granules = pd.concat([self._granules, data]).drop_duplicates(subset=["date_UTC", "tile"])

    def download_granules(self, granules: List[earthaccess.search.DataGranule]) -> List[str]:
        # Check if any of the granules have already been downloaded, and if so record the file path for that granule.
        #  Save the granules that haven't been downloaded to download them later.
        granules_to_download = []
        output_paths = []
        for granule in granules:
            date_UTC = get_date(granule["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"])
            output_file_path = join(
                self.download_directory,
                "VNP09GA",
                f"{date_UTC:%Y.%m.%d}",
                Path(granule.data_links()[0]).name
            )
            if Path(output_file_path).exists():
                output_paths.append(output_file_path)
            else:
                granules_to_download.append(granule)

        # Early exit
        if len(granules_to_download) == 0:
            logger.info("All VIIRS granules have already been downloaded")
            return output_paths

        # Make sure to remove this before we return, so we use try..finally to avoid exceptions causing issues
        temporary_download_directory = tempfile.mkdtemp(dir=join(self.download_directory, "tmp"))
        try:
            last_download_exception = None
            for _ in range(0, RETRIES):
                download_exception = None
                downloaded_granules = []

                file_paths = earthaccess.download(granules_to_download, local_path=temporary_download_directory)

                for (granule, download_file_path) in zip(granules_to_download, file_paths):
                    if isinstance(download_file_path, Exception):
                        if download_exception is None:
                            download_exception = download_file_path
                        continue
                    date_UTC = get_date(granule["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"])

                    download_file_path = Path(download_file_path)
                    output_file_path = join(
                        self.download_directory,
                        "VNP09GA",
                        f"{date_UTC:%Y.%m.%d}",
                        download_file_path.name
                    )
                    Path(output_file_path).parent.mkdir(parents=True, exist_ok=True)
                    download_file_path.rename(output_file_path)

                    output_paths.append(output_file_path)
                    downloaded_granules.append(granule)

                if download_exception is not None:
                    last_download_exception = download_exception
                    for granule in downloaded_granules:
                        granules_to_download.remove(granule)
                    logger.warning("Encountered an exception while downloading VIIRS files:", exc_info=download_exception)
                    logger.info(f"Retrying the VIIRS download with the remaining {len(granules_to_download)} granules.")
                else:
                    granules_to_download = []
                    break

            if len(granules_to_download) > 0:
                raise DownloadFailed("Error when downloading VIIRS files") from last_download_exception
        finally:
            Path(temporary_download_directory).rmdir()

        return output_paths

    def prefetch_VNP09GA(
            self,
            start_date: Union[date, str],
            end_date: Union[date, str],
            geometry: Point or Polygon or RasterGeometry = None):
        # Fetch list of granules to download
        granules = VIIRS_CMR_query(
            start_date,
            end_date,
            geometry,
        )

        self.add_granules(granules)

        self.download_granules(granules)

    def search(
            self,
            date_UTC: date,
            tile: str) -> Union[earthaccess.search.DataGranule, None]:
        subset = self._granules[(self._granules.date_UTC == date_UTC) & (self._granules.tile == tile)]
        if len(subset) > 0:
            return subset.iloc[0].granule

        granules = VIIRS_CMR_query(
            start_date=date_UTC,
            end_date=date_UTC,
            tile=tile,
        )

        if len(granules) == 0:
            return None

        if len(granules) > 0:
            logger.warning("Found more VIIRS granules than expected")

        self.add_granules(granules)

        return granules[0]

    def granule(
            self,
            date_UTC: date,
            tile: str) -> VNP09GAGranule:
        if isinstance(date_UTC, str):
            date_UTC = parser.parse(date_UTC).date()

        logger.info(f"searching VNP09GA tile {tile} date {date_UTC}")
        granule = self.search(
            date_UTC=date_UTC,
            tile=tile
        )

        if granule is None:
            raise VIIRSUnavailableError(f"VNP09GA URL not available at tile {tile} on date {date_UTC}")

        output_path = self.download_granules([granule])[0]

        output_granule = VNP09GAGranule(
            filename=output_path,
            products_directory=self.products_directory
        )

        return output_granule
