import logging
import posixpath
from os import makedirs, system, remove
from os.path import expanduser, join, exists, splitext
from shutil import move
from zipfile import ZipFile
import rasters as rt
from rasters import RasterGeometry
from timer import Timer

DEFAULT_WORKING_DIRECTORY = "."
DEFAULT_DOWNLOAD_DIRECTORY = "NLCD_download"
DEFAULT_NLCD_URL = "https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2019_land_cover_l48_20210604.zip"

logger = logging.getLogger(__name__)

class NLCDUnreachable(ConnectionError):
    pass

class NLCD:
    logger = logging.getLogger(__name__)

    def __init__(self, working_directory: str = None, download_directory: str = None, URL: str = None):
        if working_directory is None:
            working_directory = DEFAULT_WORKING_DIRECTORY

        if working_directory.startswith("~"):
            working_directory = expanduser(working_directory)

        if download_directory is None:
            download_directory = join(working_directory, DEFAULT_DOWNLOAD_DIRECTORY)

        if download_directory.startswith("~"):
            download_directory = expanduser(download_directory)

        if URL is None:
            URL = DEFAULT_NLCD_URL

        self.working_directory = working_directory
        self.download_directory = download_directory
        self.URL = URL

    def __repr__(self):
        return f'NLCD(\nworking_directory="{self.working_directory}"\ndownload_directory="{self.download_directory}"\nURL={self.URL}\n)\n'

    @property
    def zip_filename(self):
        return join(self.download_directory, posixpath.basename(self.URL))

    @property
    def granule_ID(self):
        return splitext(posixpath.basename(self.URL))[0]

    @property
    def unzipped_directory(self):
        return join(self.download_directory, self.granule_ID)

    @property
    def image_filename(self):
        return join(self.unzipped_directory, f"{self.granule_ID}.img")

    def download(self) -> str:
        if exists(self.image_filename):
            logger.info(f"NLCD file already exists: {self.image_filename}")
            return self.image_filename

        filename = self.zip_filename
        makedirs(self.download_directory, exist_ok=True)
        partial_filename = f"{filename}.download"
        command = f'wget -c -O "{partial_filename}" "{self.URL}"'
        timer = Timer()
        system(command)

        if not exists(partial_filename):
            raise NLCDUnreachable(f"unable to download URL: {self.URL}")

        move(partial_filename, filename)
        self.logger.info(f"completed download ({timer} s): {filename}")

        self.logger.info(f"unzipping {filename} -> {self.unzipped_directory}")

        with ZipFile(filename, "r") as file:
            file.extractall(self.unzipped_directory)

        if not exists(self.image_filename):
            raise NLCDUnreachable(f"failed to produce NLCD file: {self.image_filename}")

        self.logger.info(f"removing zip file: {filename}")
        remove(filename)

        return self.image_filename

    def NLCD(self, geometry: RasterGeometry):
        filename = self.download()
        image = rt.Raster.open(filename, geometry=geometry)

        return image
