""" Test of model_gen.core

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-08-13
:Copyright: 2019, Karr Lab
:License: MIT
"""

from rand_wc_model_gen.model_gen import core
import unittest

class RandomModelGeneratorTestCase(unittest.TestCase):
    def test_init(self):
        random_model = core.RandomModelGenerator()
