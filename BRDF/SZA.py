from typing import Union

import numpy as np
from rasters import Raster, RasterGeometry


def day_angle_rad_from_doy(doy):
    """
    This function calculates day angle in radians from day of year between 1 and 365.
    """
    return (2 * np.pi * (doy - 1)) / 365


def solar_dec_deg_from_day_angle_rad(day_angle_rad):
    """
    This function calculates solar declination in degrees from day angle in radians.
    """
    return (0.006918 - 0.399912 * np.cos(day_angle_rad) + 0.070257 * np.sin(day_angle_rad) - 0.006758 * np.cos(
        2 * day_angle_rad) + 0.000907 * np.sin(2 * day_angle_rad) - 0.002697 * np.cos(
        3 * day_angle_rad) + 0.00148 * np.sin(
        3 * day_angle_rad)) * (180 / np.pi)


def SZA_deg_from_lat_dec_hour(latitude, solar_dec_deg, hour):
    """
    This function calculates solar zenith angle from longitude, solar declination, and solar time.
    SZA calculated by this function matching SZA provided by MOD07 to within 0.4 degrees.

    :param latitude: latitude in degrees
    :param solar_dec_deg: solar declination in degrees
    :param hour: solar time in hours
    :return: solar zenith angle in degrees
    """
    # convert angles to radians
    latitude_rad = np.radians(latitude)
    solar_dec_rad = np.radians(solar_dec_deg)
    hour_angle_deg = hour * 15.0 - 180.0
    hour_angle_rad = np.radians(hour_angle_deg)
    sza_rad = np.arccos(
        np.sin(latitude_rad) * np.sin(solar_dec_rad) + np.cos(latitude_rad) * np.cos(solar_dec_rad) * np.cos(
            hour_angle_rad))
    sza_deg = np.degrees(sza_rad)

    return sza_deg


def calculate_SZA(day_of_year: Union[Raster, float], hour_of_day: Union[Raster, float], geometry: RasterGeometry) -> Raster:
    """
    Calculates solar zenith angle at latitude and solar apparent time.
    :param lat: latitude in degrees
    :param dec: solar declination in degrees
    :param hour: hour of day
    :return: solar zenith angle in degrees
    """
    day_angle = (2 * np.pi * (day_of_year - 1)) / 365
    lat = np.radians(geometry.lat)
    dec = np.radians((0.006918 - 0.399912 * np.cos(day_angle) + 0.070257 * np.sin(day_angle) - 0.006758 * np.cos(
        2 * day_angle) + 0.000907 * np.sin(2 * day_angle) - 0.002697 * np.cos(3 * day_angle) + 0.00148 * np.sin(
        3 * day_angle)) * (180 / np.pi))
    hour_angle = np.radians(hour_of_day * 15.0 - 180.0)
    SZA = np.degrees(np.arccos(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(hour_angle)))
    SZA = Raster(SZA, geometry=geometry)

    return SZA

