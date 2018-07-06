""" Tests API

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-03-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

import rand_wc_model_gen
import wc_model_gen
import types
import unittest


class ApiTestCase(unittest.TestCase):
    def test(self):
        self.assertIsInstance(rand_wc_model_gen, types.ModuleType)

        self.assertIsInstance(rand_wc_model_gen.config, types.ModuleType)
        self.assertIsInstance(rand_wc_model_gen.config.get_config, types.FunctionType)

        self.assertIsInstance(rand_wc_model_gen.kb_gen, types.ModuleType)
        self.assertIsInstance(rand_wc_model_gen.kb_gen.KbGenerator, type)

        self.assertIsInstance(wc_model_gen.rand_gen, types.ModuleType)
        self.assertIsInstance(wc_model_gen.rand_gen.RandomModelGenerator, type)
