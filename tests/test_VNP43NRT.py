"""
This module contains the unit tests for the VNP43NRT package.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import os
import unittest

__author__ = 'Gregory Halverson'

directory = os.path.abspath(os.path.dirname(__file__))


class TestVNP43NRT(unittest.TestCase):
    def test_import_VNP43NRT(self):
        print('testing VNP43NRT import')
        from VNP43NRT import VNP43NRT



