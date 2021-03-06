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
import tempfile
import unittest

@unittest.skip("broken_legacy")
class Test(unittest.TestCase):
    def setUp(self):
        fid, self.extra_path = tempfile.mkstemp()
        os.close(fid)

    def tearDown(self):
        os.remove(self.extra_path)

    def test_get_config(self):
        config = rand_wc_model_gen.config.get_config()
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['component']['GenomeGenerator']['num_chromosomes'], 1)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['component']['GenomeGenerator']['mean_gc_frac'], 0.5)

    def test_get_config_extra_path(self):
        with open(self.extra_path, 'w') as file:
            file.write('[rand_wc_model_gen]\n')
            file.write('    [[kb_gen]]\n')
            file.write('        [[[component]]]\n')
            file.write('            [[[[GenomeGenerator]]]]\n')
            file.write('                num_chromosomes = 5\n')

        config = rand_wc_model_gen.config.get_config(extra_path=self.extra_path)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['component']['GenomeGenerator']['num_chromosomes'], 5)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['component']['GenomeGenerator']['mean_gc_frac'], 0.5)

    def test_get_config_extra_vals(self):
        config = rand_wc_model_gen.config.get_config(extra_vals={
            'rand_wc_model_gen': {
                'kb_gen': {
                    'component': {
                        'GenomeGenerator': {
                            'num_chromosomes': 10,
                        },
                    },
                }
            }
        })
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['component']['GenomeGenerator']['num_chromosomes'], 10)
        self.assertEqual(config['rand_wc_model_gen']['kb_gen']['component']['GenomeGenerator']['mean_gc_frac'], 0.5)
