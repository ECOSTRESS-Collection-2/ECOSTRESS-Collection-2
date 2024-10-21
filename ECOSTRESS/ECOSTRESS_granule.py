import json
import logging
import os
from abc import abstractmethod
from datetime import datetime, date
from os import makedirs
from os.path import splitext, basename, abspath, expanduser, dirname, exists
from typing import Union, List

import PIL
from dateutil import parser
from matplotlib.colors import Colormap
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

import colored_logging as cl
from rasters import RasterGeometry, BBox, Point, Raster
from timer import Timer

DEFAULT_JSON_INDENT = 2

logger = logging.getLogger(__name__)


class ECOSTRESSGranule:
    _LEVEL = None
    _PRIMARY_VARIABLE = None
    _DEFAULT_UTM_CELL_SIZE = 70
    _DEFAULT_GEOGRAPHIC_CELL_SIZE = 0.0006
    _STANDARD_METADATA_GROUP = "StandardMetadata"
    _PRODUCT_METADATA_GROUP = "ProductMetadata"

    _GRANULE_PREVIEW_CMAP = "jet"

    _THUMB_SMALL_DIM = (240, 240)
    _THUMB_SMALL_QUALITY = 75
    _THUMB_LARGE_QUALITY = 30

    _GRANULE_PREVIEW_QUALITY = 60
    _GRANULE_PREVIEW_SHAPE = (1080, 1080)

    _LAYER_PREVIEW_QUALITY = 20

    def __init__(
            self,
            product_filename: str,
            *args,
            granule_preview_cmap: Union[Colormap, str] = None,
            granule_preview_shape: (int, int) = None,
            granule_preview_quality: int = None,
            product_metadata_group_name: str = None,
            standard_metadata_group_name: str = None,
            **kwargs):
        self.product_filename = product_filename
        self._geometry = None
        self._granule_preview_cmap = granule_preview_cmap
        self._granule_preview_shape = granule_preview_shape
        self._granule_preview_quality = granule_preview_quality

        if product_metadata_group_name is None:
            product_metadata_group_name = self._PRODUCT_METADATA_GROUP

        self._product_metadata_group_name = product_metadata_group_name

        if standard_metadata_group_name is None:
            standard_metadata_group_name = self._STANDARD_METADATA_GROUP

        self._standard_metadata_group_name = standard_metadata_group_name

    @property
    def primary_variable(self):
        if self._PRIMARY_VARIABLE is None:
            raise ValueError("primary variable undefined")

        return self._PRIMARY_VARIABLE

    @property
    def storage(self):
        return os.stat(self.product_filename).st_size

    @property
    def geometry(self) -> RasterGeometry:
        if self._geometry is not None:
            return self._geometry
        else:
            raise ValueError("geometry not defined")

    @property
    def bbox(self) -> BBox:
        return self.geometry.bbox.latlon

    @property
    def centroid(self) -> Point:
        return self.geometry.centroid_latlon

    @property
    def granule_name(self) -> str:
        return splitext(basename(self.product_filename))[0]

    @property
    def level(self):
        return self.granule_name.split("_")[1]

    @property
    def orbit(self) -> int:
        return int(self.granule_name.split('_')[-5])

    @property
    def scene(self) -> int:
        return int(self.granule_name.split('_')[-4])

    @property
    def product(self) -> str:
        return "_".join(self.granule_name.split("_")[1:-5])

    @property
    def timestamp(self) -> str:
        return self.granule_name.split('_')[-3]

    @property
    def build(self) -> str:
        return self.granule_name.split('_')[-2]

    @property
    def process_count(self) -> int:
        return int(self.granule_name.split('_')[-1])

    @property
    def time_UTC(self) -> datetime:
        return parser.parse(self.timestamp)

    @property
    def date_UTC(self) -> date:
        return self.time_UTC.date()

    @property
    @abstractmethod
    def boundary(self) -> Polygon:
        pass

    @property
    def boundary_WKT(self) -> str:
        return self.boundary.wkt

    @property
    @abstractmethod
    def variables(self):
        pass

    @abstractmethod
    def variable(self, variable: str) -> Raster:
        pass

    @abstractmethod
    def read_metadata_dataset(self, dataset_name: str) -> dict:
        pass

    @property
    def standard_metadata_group_name(self):
        standard_metadata_group_name = self._standard_metadata_group_name

        if standard_metadata_group_name is None:
            raise ValueError("standard metadata group not given")

        return standard_metadata_group_name

    @property
    @abstractmethod
    def standard_metadata(self) -> dict:
        pass

    @property
    @abstractmethod
    def product_metadata(self) -> dict:
        pass

    @property
    def band_specification(self) -> List[float]:
        return self.product_metadata["BandSpecification"]

    @property
    def available_bands(self) -> List[int]:
        return [band for band, wavelength in enumerate(self.band_specification) if wavelength != 0]

    @property
    def product_metadata_group_name(self):
        product_metadata_group_name = self._product_metadata_group_name

        if product_metadata_group_name is None:
            raise ValueError("product metadata group not given")

        return product_metadata_group_name

    @property
    @abstractmethod
    def product_metadata_group_path(self) -> str:
        # return f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.product_metadata_group_name}"
        pass

    @property
    @abstractmethod
    def standard_metadata_group_path(self) -> str:
        # return f"HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/{self.standard_metadata_group_name}"
        pass

    @property
    def metadata_dict(self):
        return {
            self.standard_metadata_group_name: self.standard_metadata,
            self.product_metadata_group_name: self.product_metadata
        }

    def get_metadata_JSON(self, indent: int = DEFAULT_JSON_INDENT, **kwargs):
        return json.dumps(self.metadata_dict, indent=indent, **kwargs)

    metadata_JSON = property(get_metadata_JSON)

    def write_metadata_JSON(self, filename, indent: int = DEFAULT_JSON_INDENT, **kwargs):
        filename = abspath(expanduser(filename))

        makedirs(dirname(filename), exist_ok=True)

        with open(filename, "w") as file:
            json.dump(self.metadata_dict, file, indent=indent, **kwargs)

    @property
    def granule_preview_cmap(self):
        if self._granule_preview_cmap is not None:
            return self._granule_preview_cmap

        self._granule_preview_cmap = self._GRANULE_PREVIEW_CMAP

        if self._granule_preview_cmap is None:
            raise ValueError("granule preview cmap not given")

        return self._granule_preview_cmap

    @property
    def granule_preview_shape(self):
        if self._granule_preview_shape is not None:
            return self._granule_preview_shape

        self._granule_preview_shape = self._GRANULE_PREVIEW_SHAPE

        if self._granule_preview_shape is None:
            raise ValueError("granule preview shape not given")

        return self._granule_preview_shape

    @property
    def granule_preview_quality(self):
        if self._granule_preview_quality is not None:
            return self._granule_preview_quality

        self._granule_preview_quality = self._GRANULE_PREVIEW_QUALITY

        if self._granule_preview_quality is None:
            raise ValueError("granule preview quality not given")

        return self._granule_preview_quality

    def get_browse_image(
            self,
            cmap: Union[Colormap, str] = None,
            mode: str = "RGB") -> PIL.Image.Image:

        if cmap is None:
            cmap = self.granule_preview_cmap

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        image = self.variable(self.primary_variable)

        browse_image = image.percentilecut.resize(self.granule_preview_shape, resampling="nearest").to_pillow(
            cmap=cmap,
            mode=mode
        )

        return browse_image

    browse_image = property(get_browse_image)

    def write_browse_image(self, PNG_filename: str, cmap: Union[Colormap, str] = None):
        if cmap is None:
            cmap = self.granule_preview_cmap

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        PNG_filename = abspath(expanduser(PNG_filename))
        makedirs(dirname(PNG_filename), exist_ok=True)
        shape = self.granule_preview_shape
        quality = self.granule_preview_quality
        logger.info(
            f"started writing PNG browse image with cmap {cl.name(cmap.name)} shape {cl.val(shape)} and quality {cl.val(quality)}: {cl.file(PNG_filename)}")
        timer = Timer()

        # image = self.variable(self.primary_variable)
        # print("self.browse_cmap")
        # print(self.browse_cmap)
        # browse_image = image.percentilecut.resize(self.granule_preview_shape, resampling="nearest").to_pillow(
        #     cmap=cmap)

        browse_image = self.get_browse_image(cmap=cmap)
        browse_image.save(PNG_filename, format="png", quality=self.granule_preview_quality)
        logger.info(f"finished writing PNG browse image ({cl.time(timer)})")

        if not exists(PNG_filename):
            raise IOError(f"unable to create PNG browse image: {PNG_filename}")
