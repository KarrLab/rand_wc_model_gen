""" Test of model_gen.core

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-08-13
:Copyright: 2019, Karr Lab
:License: MIT
"""

from rand_wc_model_gen.model_gen import core
import unittest

class RandModelGenTestCase(unittest.TestCase):
    def test_init(self):
        random_model_generator_1 = core.RandModelGen()
        self.assertEqual(random_model_generator_1.options, {'id':None, 'name':None, 'version':None})

        random_model_generator_2 = core.RandModelGen(options={'id':'test_rand', 'name':'test random model', 'version':'0.0'})
        self.assertEqual(random_model_generator_2.options, {'id':'test_rand', 'name':'test random model', 'version':'0.0'})

    def test_run(self):
        random_model_1 = core.RandModelGen().run()
        self.assertIsNone(random_model_1.id)
        self.assertIsNone(random_model_1.name)
        self.assertIsNone(random_model_1.version)

        random_model_2 = core.RandModelGen(options={'id':'test_rand', 'name':'test random model', 'version':'0.0'}).run()
        self.assertEqual(random_model_2.id, 'test_rand')
        self.assertEqual(random_model_2.name, 'test random model')
        self.assertEqual(random_model_2.version, '0.0')
