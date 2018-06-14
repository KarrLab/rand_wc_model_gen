""" Tests of transcription submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import kb_gen
from rand_wc_model_gen.model_gen import transcription
import numpy
import scipy
import unittest
import wc_kb
import wc_lang


class TranscriptionSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = kb_gen.KbGenerator(options={
            'component': {
                'PropertiesGenerator': {
                    'mean_volume': 1e-15,
                    'mean_doubling_time': 1000.,
                },
                'ChromosomesGenesTusGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 100.,
                    'mean_gene_len': 100.,
                },
                'MetabolitesGenerator': {
                },
                'RnaGenerator': {
                    'mean_copy_number': 10.,
                    'mean_half_life': 120.,
                },
            },
        }).run()
        cell = kb.cell
        rna = cell.species_types.get(__type=wc_kb.RnaSpeciesType)

        model = wc_lang.Model()
        gen = transcription.TranscriptionSubmodelGenerator(kb, model, options={})
        gen.run()

        submodel = model.submodels.get_one(id='transcription')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

        # check species types and species generated
        atp = model.species_types.get_one(id='atp')
        atp_cytosol = atp.species.get_one(compartment=cytosol)
        self.assertEqual(atp_cytosol.concentration.units, 'M')

        concs = []
        for species_type in model.species_types:
            if species_type.id.startswith('rna_'):
                species = species_type.species.get_one(compartment=cytosol)
                concs.append(species.concentration.value)
                self.assertEqual(species.concentration.units, 'M')
        numpy.testing.assert_almost_equal(numpy.mean(concs), 10. / scipy.constants.Avogadro / 1e-15, decimal=2)

        # check reactions generated
        self.assertEqual(len(submodel.reactions), len(rna))
        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=atp).coefficient
            + submodel.reactions[0].participants.get_one(species=ctp).coefficient
            + submodel.reactions[0].participants.get_one(species=gtp).coefficient
            + submodel.reactions[0].participants.get_one(species=utp).coefficient,
            -rna[0].get_len())
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=ppi).coefficient,
            rna[0].get_len())
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h2o).coefficient,
            rna[0].get_len()-1)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h).coefficient,
            -(rna[0].get_len()-1))

        # check rate laws
        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertEqual(rl.direction.name, 'forward')
            self.assertEqual(rl.equation.expression, 'k_cat')
            self.assertEqual(rl.equation.modifiers, [])
            self.assertEqual(rl.equation.parameters, [])            
            self.assertEqual(rl.k_m, None)

        k_cats = [rxn.rate_laws[0].k_cat for rxn in submodel.reactions]
        numpy.testing.assert_almost_equal(numpy.mean(k_cats), 10. * numpy.log(2.) * (1. / 120.), decimal=2)
