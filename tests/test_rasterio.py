"""
This module contains the unit tests for the rasterio package.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import os
import unittest

__author__ = 'Gregory Halverson'

directory = os.path.abspath(os.path.dirname(__file__))


class TestRasterio(unittest.TestCase):
    def test_import_rasterio(self):
        print('testing rasterio import')
        import rasterio as rio
        print("rasterio version: {}".format(rio.__version__))


