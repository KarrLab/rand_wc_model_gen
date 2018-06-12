""" Tests of model generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import kb_gen
from rand_wc_model_gen import model_gen
import unittest
import wc_lang


class ModelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = kb_gen.KbGenerator(options={
            'component': {
                'ChromosomesGenesGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 100,
                    'mean_gene_len': 100,
                },
            },
        }).run()
        model = model_gen.ModelGenerator(kb).run()

        self.assertIsInstance(model.submodels.get_one(id='transcription'), wc_lang.Submodel)
