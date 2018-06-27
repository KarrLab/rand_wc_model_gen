""" Tests of knowledge base generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import kb_gen
import obj_model
import unittest
import wc_kb
import wc_utils.util.string


class KbGeneratorTestCase(unittest.TestCase):
    def test(self):
        gen = kb_gen.KbGenerator(options={
            'component': {
                'PropertiesGenerator': {
                    'mean_volume': 1e-15,
                    'mean_doubling_time': 1000.,
                },
                'GenomeGenerator': {
                    'num_chromosomes': 10,
                    'mean_num_genes': 100,

                    'mean_num_sRNA': 10,
                    'mean_num_rRNA': 10,
                    'mean_num_tRNA': 10}}})
        
        kb = gen.run()
        self.assertEqual(len(kb.cell.species_types.get(
            __type=wc_kb.DnaSpeciesType)), 10)

        errors = obj_model.Validator().run(kb, get_related=True)
        self.assertEqual(
            errors, None, msg=wc_utils.util.string.indent_forest(errors))

    def test_clean_and_validate_options(self):
        gen = kb_gen.KbGenerator()

        gen.options = {'seed': None}
        gen.clean_and_validate_options()

        gen.options = {'seed': int(2)}
        gen.clean_and_validate_options()

        gen.options = {'seed': 'ABC'}
        with self.assertRaises(ValueError):
            gen.clean_and_validate_options()
