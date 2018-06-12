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
                'ChromosomesGenesGenerator': {
                    'num_chromosomes': 10,
                }
            }
        })
        kb = gen.run()
        self.assertEqual(len(kb.cell.species_types.get(__type=wc_kb.DnaSpeciesType)), 10)
