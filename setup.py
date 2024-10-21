from setuptools import setup, find_packages

try:
    import pybind11_cmake
except ImportError:
    print("pybind11-cmake must be installed."
          "Try \n \t pip install pybind11_cmake")
    import sys
    sys.exit()

from pybind11_cmake import CMakeExtension, CMakeBuild

from os.path import join, abspath, dirname

__author__ = "Gregory Halverson"
AUTHOR_EMAIL = "Gregory.H.Halverson@jpl.nasa.gov"
# PACKAGE_NAME = "ECOSTRESS"
# DESCRIPTION = "Spatially and Temporally Adaptive Remote Sensing (STARS) Data Fusion Model"

def version():
    with open(join(abspath(dirname(__file__)), "ECOSTRESS", "PGEVersion.txt"), "r") as file:
        return file.read()

setup(
    # name=PACKAGE_NAME,
    version=version(),
    author=__author__,
    author_email=AUTHOR_EMAIL,
    # description=DESCRIPTION,
    long_description='',
    setup_requires=['pybind11_cmake'],
    ext_modules=[CMakeExtension('CPPSTARS')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    packages=find_packages() + ["STARS_jl", "STARS_jl.src", "VNP43NRT_jl", "VNP43NRT_jl.src"],
    package_data={'': ["*", "STARS_jl/*", "STARS_jl/src/*", "VNP43NRT_jl/*", "VNP43NRT_jl/src/*"]}
)
