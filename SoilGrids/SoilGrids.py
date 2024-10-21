import logging
import posixpath
from os import makedirs, system
from os.path import abspath, expanduser, exists, dirname, join
from shutil import move
from time import perf_counter
import rasters as rt
from rasters import RasterGeometry, Raster
import numpy as np

DEFAULT_RESAMPLING = "cubic"
DEFAULT_WORKING_DIRECTORY = "."
DEFAULT_DOWNLOAD_DIRECTORY = "SoilGrids_download"

class SoilGrids:
    # FC_URL = "https://zenodo.org/record/2784001/files/sol_watercontent.33kPa_usda.4b1c_m_250m_b100..100cm_1950..2017_v0.1.tif"
    FC_URL = "https://zenodo.org/record/2784001/files/sol_watercontent.33kPa_usda.4b1c_m_250m_b0..0cm_1950..2017_v0.1.tif"
    # WP_URL = "https://zenodo.org/record/2784001/files/sol_watercontent.1500kPa_usda.3c2a1a_m_250m_b100..100cm_1950..2017_v0.1.tif"
    WP_URL = "https://zenodo.org/record/2784001/files/sol_watercontent.1500kPa_usda.3c2a1a_m_250m_b0..0cm_1950..2017_v0.1.tif"

    logger = logging.getLogger(__name__)

    def __init__(self, working_directory: str = None, source_directory: str = None, resampling: str = DEFAULT_RESAMPLING):
        if working_directory is None:
            working_directory = abspath(expanduser(DEFAULT_WORKING_DIRECTORY))

        if source_directory is None:
            source_directory = join(working_directory, DEFAULT_DOWNLOAD_DIRECTORY)

        self.working_directory = abspath(expanduser(working_directory))
        self.source_directory = abspath(expanduser(source_directory))
        self.resampling = resampling

    def __repr__(self) -> str:
        return f'SoilGrids(source_directory="{self.source_directory}")'

    @property
    def FC_filename(self):
        return join(self.source_directory, posixpath.basename(self.FC_URL))

    @property
    def WP_filename(self):
        return join(self.source_directory, posixpath.basename(self.WP_URL))

    def download_file(self, URL: str, filename: str) -> str:
        if exists(filename):
            self.logger.info(f"file already downloaded: {filename}")
            return filename

        self.logger.info(f"downloading: {URL} -> {filename}")
        directory = dirname(filename)
        makedirs(directory, exist_ok=True)
        partial_filename = f"{filename}.download"
        command = f'wget -c -O "{partial_filename}" "{URL}"'
        download_start = perf_counter()
        system(command)
        download_end = perf_counter()
        download_duration = download_end - download_start
        self.logger.info(f"completed download in {download_duration:0.2f} seconds: {filename}")

        if not exists(partial_filename):
            raise IOError(f"unable to download URL: {URL}")

        move(partial_filename, filename)

        return filename

    def FC(self,
           geometry: RasterGeometry = None,
           resampling: str = None) -> Raster:
        if resampling is None:
            resampling = self.resampling

        URL = self.FC_URL
        filename = self.FC_filename

        if not exists(filename):
            self.download_file(URL, filename)

        image = rt.Raster.open(
            filename=filename,
            geometry=geometry,
            resampling=resampling
        )

        image = rt.where(image == 255, np.nan, image)
        image.nodata = np.nan
        image = rt.clip(image / 100, 0, 1)

        return image

    def WP(self,
           geometry: RasterGeometry = None,
           resampling: str = None) -> Raster:
        if resampling is None:
            resampling = self.resampling

        URL = self.WP_URL
        filename = self.WP_filename

        if not exists(filename):
            self.download_file(URL, filename)

        image = rt.Raster.open(
            filename=filename,
            geometry=geometry,
            resampling=resampling
        )

        image = rt.where(image == 255, np.nan, image)
        image.nodata = np.nan
        image = rt.clip(image / 100, 0, 1)

        return image
