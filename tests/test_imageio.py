"""
This module contains the unit tests for imageio compatibility.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import numpy as np
import os
import unittest

__author__ = 'Gregory Halverson'


class TestImageIO(unittest.TestCase):
    def test_import_imageio(self):
        import imageio
