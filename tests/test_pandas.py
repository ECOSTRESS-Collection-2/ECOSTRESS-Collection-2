"""
This module contains the unit tests for the pandas package.

Developed by Gregory Halverson in the Jet Propulsion Laboratory Year-Round Internship Program (Columbus Technologies and Services), in coordination with the ECOSTRESS mission and master's thesis studies at California State University, Northridge.
"""

import unittest
from os.path import join, abspath, dirname

__author__ = 'Gregory Halverson'


class TestPandas(unittest.TestCase):
    def test_import_pandas(self):
        print('testing pandas import')
        import pandas as pd
        print("pandas version: {}".format(pd.__version__))

    def test_sort_values(self):
        print('testing pandas.sort_values')
        import pandas as pd

        region_queue = pd.read_csv(join(abspath(dirname(__file__)), 'region_queue.csv'))
        region_queue.sort_values(
            ['target_date', 'target_region'],
            ascending=[True, True],
            inplace=True
        )
