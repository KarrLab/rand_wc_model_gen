""" Test model generation

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-01-08
:Copyright: 2018, Karr Lab
:License: MIT
"""
import unittest
import os

from wc_lang.io import Reader
#from wc_utils.workbook.io import read as read_workbook
from random_wc_model_generator.core import GenerateModel, main


class TestGenerateModel(unittest.TestCase):

    def setUp(self):
        fixtures = os.path.join(os.path.dirname(__file__), 'fixtures')
        self.CONIG_FILENAME = os.path.join(fixtures, 'config.yml')
        self.MODEL_FILENAME = os.path.join(fixtures, 'wc_model')

    def generate_model(self, config_file, output_file, model_format=None):
        args = [config_file, output_file]
        if model_format is not None:
            args.append("--model_format {}".format(model_format))
        return main(args)

    def test_generate_model(self):
        generate_model = self.generate_model(self.CONIG_FILENAME, self.MODEL_FILENAME + '.xlsx')
        self.assertEqual(generate_model.config['genes'], 100)

    def test_metabolism_core(self):
        metabolism_core = os.path.join(os.path.dirname(__file__),
            '../random_wc_model_generator/data/fixtures', 'metabolism.xlsx')
        model = Reader().run(metabolism_core)
