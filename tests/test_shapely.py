"""
This module contains the unit tests for the shapely package.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import os
import unittest

__author__ = 'Gregory Halverson'

directory = os.path.abspath(os.path.dirname(__file__))


class TestShapely(unittest.TestCase):
    def test_import_shapely(self):
        print('testing shapely import')
        import shapely
        print("shapely version: {}".format(shapely.__version__))

    def test_polygon_creation(self):
        from sentinel_tile_grid import sentinel_tile_grid
        grid = sentinel_tile_grid.grid("13SFU")
