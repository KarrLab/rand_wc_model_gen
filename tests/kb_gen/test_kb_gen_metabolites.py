""" Tests of generation of chromosomes and genes

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen.kb_gen import metabolites
import numpy
import unittest
import wc_kb


class MetabolitesGeneratorTestCase(unittest.TestCase):
    def test_run(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        gen = metabolites.MetabolitesGenerator(kb, options={
            'mean_ntp_conc': 10.,
            'mean_h2o_conc': 20.,
            'mean_ph': 7.,
        })
        gen.run()

        atp = cell.species_types.get_one(__type=wc_kb.MetaboliteSpeciesType, id='atp')
        ppi = cell.species_types.get_one(__type=wc_kb.MetaboliteSpeciesType, id='ppi')
        h2o = cell.species_types.get_one(__type=wc_kb.MetaboliteSpeciesType, id='h2o')
        h = cell.species_types.get_one(__type=wc_kb.MetaboliteSpeciesType, id='h')
        self.assertEqual(atp.concentration, 10.)
        self.assertEqual(ppi.concentration, 10.)
        self.assertEqual(h2o.concentration, 20.)
        self.assertEqual(h.concentration, numpy.power(10., -7.))
