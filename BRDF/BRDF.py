from os.path import join, abspath, dirname
import warnings
import numpy as np
from numpy import where, isnan, logical_not, nanmean, nan, logical_or, loadtxt, digitize, arange


def statistical_radiative_transport(SZA, AOT):
    filename = join(abspath(dirname(__file__)), 'statistical_radiative_transport.txt')

    data = loadtxt(filename)

    sza_index = digitize(SZA, arange(0, 90, 1)[:-1])
    aot_index = digitize(AOT, arange(0, 1, 0.02)[:-1])

    SRT = data[sza_index, aot_index]

    return SRT


def bidirectional_reflectance(white_sky_albedo, black_sky_albedo, SZA, AOT):
    warnings.filterwarnings('ignore')

    # constrain aerosol optical thickness
    # AOT = where(AOT >= 0.98, 0.97, AOT)
    # AOT = where(AOT < 0.02, 0.02, AOT)
    # AOT = where(isnan(AOT), 0.1, AOT)
    AOT = np.clip(AOT, 0.1, 0.97)
    AOT_mean = np.nanmean(AOT)
    AOT = np.where(np.isnan(AOT), AOT_mean, AOT)

    # constrain SZA
    SZA = where(
        logical_or(SZA <= 0, SZA > 90),
        nanmean(where(
            logical_not(logical_or(SZA <= 0, SZA > 90)),
            SZA,
            nan
        )),
        SZA
    )

    warnings.resetwarnings()

    # gap-fill SZA with constant 45 degrees
    SZA = where(isnan(SZA), 45, SZA)

    # lookup statistical radiative transport
    SRT = statistical_radiative_transport(SZA, AOT)

    warnings.filterwarnings('ignore')

    # cloud filter albedo
    blue_sky_albedo = white_sky_albedo * SRT + black_sky_albedo * (1 - SRT)

    warnings.resetwarnings()

    return blue_sky_albedo
