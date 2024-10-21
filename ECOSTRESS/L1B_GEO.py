import logging
import os
from os.path import abspath, expanduser, exists

import h5py
import numpy as np
from numpy import array
from shapely.geometry import Polygon
from six import string_types

from ECOSTRESS.ECOSTRESS_HDF5_granule import ECOSTRESSHDF5Granule
from ECOSTRESS.exit_codes import InputFilesInaccessible, UncalibratedGeolocation, InvalidCoordinates
from HDFEOS5.HDF5 import HDF5
from rasters import Raster, RasterGeolocation, WGS84

logger = logging.getLogger(__name__)


class L1BGEO(ECOSTRESSHDF5Granule):
    _LEVEL = "L1B"
    _L1B_GEO_DATA_GROUP = "Geolocation"
    _DATA_GROUP = _L1B_GEO_DATA_GROUP
    _L1B_GEO_METADATA_GROUP = "L1GEOMetadata"
    _PRODUCT_METADATA_GROUP = _L1B_GEO_METADATA_GROUP
    _LAT_NAME = "latitude"
    _LON_NAME = "longitude"
    _CORRECTED_GEOLOCATION = "OrbitCorrectionPerformed"

    def __init__(self, L1B_GEO_filename: str, *args, **kwargs):
        if L1B_GEO_filename is None:
            raise ValueError("L1B GEO filename not given")

        L1B_GEO_filename = abspath(expanduser(L1B_GEO_filename))

        if not exists(L1B_GEO_filename):
            raise InputFilesInaccessible(f"L1B GEO file not found: {L1B_GEO_filename}")

        ECOSTRESSHDF5Granule.__init__(self, product_filename=L1B_GEO_filename)
        self.L1B_GEO_filename = L1B_GEO_filename
        self._geolocation = None

    def __repr__(self):
        return f'L1B_GEO("{self.L1B_GEO_filename}")'

    @property
    def land_percent(self) -> float:
        with h5py.File(self.L1B_GEO_filename, "r") as f:
            return float(f['L1GEOMetadata/OverAllLandFraction'][()])

    @property
    def cloud(self) -> None:
        return None

    @property
    def view_zenith(self) -> Raster:
        with h5py.File(self.L1B_GEO_filename, "r") as file:
            data = np.array(file["Geolocation/view_zenith"])

        image = Raster(data, geometry=self.geolocation)

        return image

    @property
    def height(self) -> Raster:
        with h5py.File(self.L1B_GEO_filename, "r") as file:
            data = np.array(file["Geolocation/height"])

        image = Raster(data, geometry=self.geolocation)

        return image

    @property
    def geolocation(self) -> RasterGeolocation:
        """
        This function reads the geolocation arrays from an ECOSTRES L1B GEO product as a CoordinateField object.
        :return: <RasterGeolocation> geolocation as CoordinateField object
        """
        if self._geolocation is not None:
            return self._geolocation

        if not isinstance(self.L1B_GEO_filename, string_types):
            raise ValueError("invalid ECOSTRESS L1B GEO filename passed when attempting to read geolocation")

        if not exists(self.L1B_GEO_filename):
            raise InputFilesInaccessible(f"L1B GEO product {self.L1B_GEO_filename} does not exist")

        if not os.access(self.L1B_GEO_filename, os.R_OK):
            raise InputFilesInaccessible(f"L1B GEO product {self.L1B_GEO_filename} does not have read permission")

        lon = None
        lat = None

        with HDF5(self.L1B_GEO_filename, 'r') as file:
            if file is None:
                raise InputFilesInaccessible(f"L1B GEO opened as null: {self.L1B_GEO_filename}")

            top_level_keys = list(file.keys())

            if self._L1B_GEO_DATA_GROUP not in top_level_keys:
                raise InputFilesInaccessible(
                    f"geolocation group {self._L1B_GEO_DATA_GROUP} "
                    f"not available in {', '.join(top_level_keys)} "
                    f"top level keys of {self.L1B_GEO_filename}"
                )

            if self._L1B_GEO_METADATA_GROUP not in top_level_keys:
                raise InputFilesInaccessible(
                    f"L1B GEO product metadata group {self._L1B_GEO_METADATA_GROUP} "
                    f"not available in {', '.join(top_level_keys)} "
                    f"top level keys of {self.L1B_GEO_filename}"
                )

            if file[self._L1B_GEO_METADATA_GROUP][self._CORRECTED_GEOLOCATION][()] == "False":
                raise UncalibratedGeolocation(f"L1B GEO product {self.L1B_GEO_filename} is uncalibrated")

            try:
                group = file.get(self._L1B_GEO_DATA_GROUP)
            except Exception as e:
                logger.error(e)
                raise InputFilesInaccessible(
                    f"group {self._L1B_GEO_DATA_GROUP} could not be read from L1B GEO file {self.L1B_GEO_filename}: {e}")

            try:
                lat = array(group.get(self._LAT_NAME))
            except Exception as e:
                logger.error(e)
                raise InputFilesInaccessible(
                    f"latitude variable {self._LAT_NAME} could not be read from L1B GEO file {self.L1B_GEO_filename} group {self._L1B_GEO_DATA_GROUP}")
            try:
                lon = array(group.get(self._LON_NAME))
            except Exception as e:
                logger.error(e)

                raise InputFilesInaccessible(
                    f"longitude variable {self._LON_NAME} "
                    f"could not be read from L1B GEO file {self.L1B_GEO_filename} "
                    f"group {self._L1B_GEO_DATA_GROUP}"
                )

        if lat is None or lon is None:
            raise InputFilesInaccessible(f"missing geolocation in L1B GEO file: {self.L1B_GEO_filename}")

        if np.any(np.logical_or(lat < -90, lat > 90)):
            raise InvalidCoordinates(f"L1B GEO product contains invalid latitude: {self.L1B_GEO_filename}")

        if np.any(np.logical_or(lon < -180, lat > 180)):
            raise InvalidCoordinates(f"L1B GEO product contains invalid latitude: {self.L1B_GEO_filename}")

        xmin = np.nanmin(lon)
        xmax = np.nanmax(lon)
        width = xmax - xmin
        ymin = np.nanmin(lat)
        ymax = np.nanmax(lat)
        height = ymax - ymin

        if width > 10 or height > 10:
            raise InvalidCoordinates(f"invalid coordinate field width: {width} height: {height}")

        self._geolocation = RasterGeolocation(lon, lat, WGS84)

        return self._geolocation

    @property
    def geometry(self) -> RasterGeolocation:
        return self.geolocation

    @property
    def corner_polygon(self) -> Polygon:
        """
        This function reads the bounding polygon of an ECOSTRESS L1B GEO scene.
        :return: <Polygon> bounding polygon as shapely Polygon object
        """
        return self.geolocation.corner_polygon
