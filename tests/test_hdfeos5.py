"""
This module contains the unit tests for HDF-EOS5.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import os
import unittest

__author__ = 'Gregory Halverson'


class TestHDFEOS5(unittest.TestCase):
    def test_import_he5py(self):
        import he5py

