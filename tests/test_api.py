""" Tests API

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-03-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

import random_wc_model_generator
import types
import unittest


class ApiTestCase(unittest.TestCase):
    def test(self):
        self.assertIsInstance(random_wc_model_generator, types.ModuleType)
        self.assertIsInstance(random_wc_model_generator.CreateWcLangModel, type)
        self.assertIsInstance(random_wc_model_generator.enrich_polymers, types.ModuleType)
