from os.path import join, abspath, dirname
from .ECOSTRESS import *

with open(join(abspath(dirname(__file__)), "version.txt")) as f:
    version = f.read()

with open(join(abspath(dirname(__file__)), 'PGEVersion.txt')) as f:
    PGEVersion = f.read().strip()

__version__ = version
__author__ = "Gregory H. Halverson"
