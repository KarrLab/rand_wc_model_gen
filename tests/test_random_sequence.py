""" Test model generation

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-01-08
:Copyright: 2018, Karr Lab
:License: MIT
"""
import unittest

import RandomSeqGen from random_sequence

class TestGenerateModel(unittest.TestCase):

    '''
    def setUp(self):
        fixtures = os.path.join(os.path.dirname(__file__), 'fixtures')
        self.CONIG_FILENAME = os.path.join(fixtures, 'config.yml')
        self.MODEL_FILENAME = os.path.join(fixtures, 'wc_model')

    def generate_model(self, config_file, output_file, model_format=None):
        args = [config_file, output_file]
        if model_format is not None:
            args.append("--model_format {}".format(model_format))
        return main(args)
    '''

    def test_random_sequence(self):
        pass

    def test_make_genes(self):
        pass

# pytest tests/test_core.py 
