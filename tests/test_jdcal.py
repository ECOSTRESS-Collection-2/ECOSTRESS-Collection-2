"""
This module contains the unit tests for the gdal package.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import os
import unittest

__author__ = 'Gregory Halverson'

directory = os.path.abspath(os.path.dirname(__file__))


class TestJDCal(unittest.TestCase):
    def test_import_jdcal(self):
        print('testing jdcal import')
        import jdcal
        print("jdcal version: {}".format(jdcal.__version__))
