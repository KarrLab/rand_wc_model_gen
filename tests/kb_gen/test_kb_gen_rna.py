""" Tests of RNA generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import kb_gen
from wc_utils.util.units import unit_registry
import numpy
import scipy
import unittest
import wc_kb

@unittest.skip("broken_legacy")
class RnaGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        cell.properties.create(id='mean_volume', value=1, units=unit_registry.parse_units('L'))
        cytosol = cell.compartments.get_one(id='c')

        tus = []
        for i_tu in range(1000):
            tu = cell.loci.create(__type=wc_kb.prokaryote.TranscriptionUnitLocus,
                                  id='tu_{}'.format(i_tu + 1),
                                  name='Transcription unit {}'.format(i_tu + 1),
                                  start=10 + 20 * (i_tu), end=20 + 20 * (i_tu), strand=wc_kb.core.PolymerStrand.positive)
            tu.genes.create(id='gene_{}'.format(i_tu + 1), type=wc_kb.core.GeneType.mRna)
            tus.append(tu)

        gen = kb_gen.rna.RnaGenerator(kb, options={
            'mean_copy_number': 10.,
            'mean_half_life': 120.,
        })
        gen.run()

        rnas = cell.species_types.get(__type=wc_kb.prokaryote.RnaSpeciesType)
        self.assertEqual(len(rnas), 1000)

        self.assertEqual(rnas[0].transcription_units, [tus[0]])
        self.assertEqual(rnas[0].name, 'RNA 1')
        self.assertEqual(rnas[0].type, wc_kb.core.RnaType.mRna)

        concs = [rna.species.get_one(compartment=cytosol).concentration.value for rna in rnas]
        half_lives = [rna.half_life for rna in rnas]
        self.assertAlmostEqual(numpy.mean(concs), 10. / scipy.constants.Avogadro / 1., delta=5. * numpy.sqrt(10. / 1000.))
        self.assertAlmostEqual(numpy.mean(half_lives), 120., delta=5. * numpy.sqrt(120. / 1000.))
