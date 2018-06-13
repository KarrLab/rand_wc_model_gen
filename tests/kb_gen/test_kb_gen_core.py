""" Tests of knowledge base generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import kb_gen
import unittest
import wc_kb


class KbGeneratorTestCase(unittest.TestCase):
    def test(self):
        gen = kb_gen.KbGenerator(options={
            'component': {
                'ChromosomesGenesTusGenerator': {
                    'num_chromosomes': 10,
                    'mean_num_genes': 100,
                }
            }
        })

        kb = gen.run()
        self.assertEqual(len(kb.cell.species_types.get(__type=wc_kb.DnaSpeciesType)), 10)

    def test_clean_and_validate_options(self):
        gen = kb_gen.KbGenerator()

        gen.options = {'seed': None}
        gen.clean_and_validate_options()

        gen.options = {'seed': int(2)}
        gen.clean_and_validate_options()

        gen.options = {'seed': 'ABC'}
        with self.assertRaises(ValueError):
            gen.clean_and_validate_options()
