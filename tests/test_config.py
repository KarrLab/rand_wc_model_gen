""" Tests of the configuration

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import rand_wc_model_gen
import rand_wc_model_gen.config
import os
import pkg_resources
import tempfile
import unittest


class Test(unittest.TestCase):
    def setUp(self):
        fid, self.extra_path = tempfile.mkstemp()
        os.close(fid)

    def tearDown(self):
        os.remove(self.extra_path)

    def test_get_config(self):
        config = rand_wc_model_gen.config.get_config()
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['genome_len'], 1000)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['num_genes'], 10)

    def test_get_config_extra_path(self):
        with open(self.extra_path, 'w') as file:
            file.write('[rand_wc_model_gen]\n')
            file.write('    [[kb_gen]]\n')
            file.write('        genome_len = 2500\n')

        config = rand_wc_model_gen.config.get_config(extra_path=self.extra_path)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['genome_len'], 2500)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['num_genes'], 10)

    def test_get_config_extra_vals(self):
        config = rand_wc_model_gen.config.get_config(extra_vals={
            'rand_wc_model_gen': {
                'kb_gen': {
                    'genome_len': 2000,
                }
            }
        })
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['genome_len'], 2000)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['num_genes'], 10)
