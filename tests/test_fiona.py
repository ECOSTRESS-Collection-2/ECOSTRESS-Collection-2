"""
This module contains the unit tests for the fiona package.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import os
import unittest

__author__ = 'Gregory Halverson'

directory = os.path.abspath(os.path.dirname(__file__))


class TestFiona(unittest.TestCase):
    def test_import_fiona(self):
        print('testing fiona import')
        import fiona
        print("fiona version: {}".format(fiona.__version__))
        print("gdal version: {}".format(fiona.__gdal_version__))
