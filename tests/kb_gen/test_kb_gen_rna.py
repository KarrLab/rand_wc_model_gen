""" Tests of RNA generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""

import numpy
from rand_wc_model_gen import kb_gen
import scipy
import unittest
import wc_kb


class RnaGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        cell.properties.create(id='mean_volume', value=1, units='L')

        for i_gene in range(1000):
            cell.loci.create(__type=wc_kb.GeneLocus,
                             id='gene_{}'.format(i_gene + 1),
                             start=10 + 20 * (i_gene), end=20 + 20 * (i_gene), strand=wc_kb.PolymerStrand.positive,
                             type=wc_kb.GeneType.mRna)

        gen = kb_gen.rna.RnaGenerator(kb, options={
            'mean_copy_number': 10,
            'mean_half_life': 120,
        })
        gen.run()

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        self.assertEqual(len(rnas), 1000)

        concs = [rna.concentration for rna in rnas]
        half_lives = [rna.half_life for rna in rnas]
        self.assertAlmostEqual(numpy.mean(concs), 10 / scipy.constants.Avogadro / 1, delta=5 * numpy.sqrt(10 / 1000))
        self.assertAlmostEqual(numpy.mean(half_lives), 120, delta=5 * numpy.sqrt(120 / 1000))
