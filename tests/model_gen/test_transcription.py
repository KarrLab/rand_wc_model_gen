""" Tests of transcription submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import kb_gen
from rand_wc_model_gen.model_gen import transcription
import unittest
import wc_kb
import wc_lang


class TranscriptionSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = kb_gen.KbGenerator(options={
            'component': {
                'ChromosomesGenesGenerator': {
                    'num_chromosomes': 1,
                    'avg_num_genes': 100,
                    'avg_gene_len': 100,
                },
            },
        }).run()
        cell = kb.cell

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
        self.assertEqual(atp_cytosol.serialize(), 'atp[c]')

        # check reactions generated
        genes = cell.loci.get(__type=wc_kb.GeneLocus)
        self.assertEqual(len(submodel.reactions), len(genes))
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
            -genes[0].get_len())
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=ppi).coefficient,
            genes[0].get_len())
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h2o).coefficient,
            genes[0].get_len()-1)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h).coefficient,
            -(genes[0].get_len()-1))
