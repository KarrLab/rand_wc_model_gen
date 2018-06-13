""" Generator for transcription submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import numpy
import scipy
import wc_kb
import wc_lang
import wc_model_gen


class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for transcription submodel """

    def gen_compartments(self):
        model = self.model
        cytosol = model.compartments.get_or_create(id='c')
        cytosol.name = 'cytosol'

    def gen_species(self):
        """ Generate species associated with submodel """
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

        # get or create metabolite species
        kb_atp = cell.species_types.get_one(id='atp')
        kb_ctp = cell.species_types.get_one(id='ctp')
        kb_gtp = cell.species_types.get_one(id='gtp')
        kb_utp = cell.species_types.get_one(id='utp')
        kb_ppi = cell.species_types.get_one(id='ppi')
        kb_h2o = cell.species_types.get_one(id='h2o')
        kb_h = cell.species_types.get_one(id='h')

        model_atp = model.species_types.get_or_create(id='atp')
        model_ctp = model.species_types.get_or_create(id='ctp')
        model_gtp = model.species_types.get_or_create(id='gtp')
        model_utp = model.species_types.get_or_create(id='utp')
        model_ppi = model.species_types.get_or_create(id='ppi')
        model_h2o = model.species_types.get_or_create(id='h2o')
        model_h = model.species_types.get_or_create(id='h')

        model_atp_c = model_atp.species.get_or_create(compartment=cytosol)
        model_ctp_c = model_ctp.species.get_or_create(compartment=cytosol)
        model_gtp_c = model_gtp.species.get_or_create(compartment=cytosol)
        model_utp_c = model_utp.species.get_or_create(compartment=cytosol)
        model_ppi_c = model_ppi.species.get_or_create(compartment=cytosol)
        model_h2o_c = model_h2o.species.get_or_create(compartment=cytosol)
        model_h_c = model_h.species.get_or_create(compartment=cytosol)

        model_atp_c.concentration = wc_lang.Concentration(value=kb_atp.concentration, units='M')
        model_ctp_c.concentration = wc_lang.Concentration(value=kb_ctp.concentration, units='M')
        model_gtp_c.concentration = wc_lang.Concentration(value=kb_gtp.concentration, units='M')
        model_utp_c.concentration = wc_lang.Concentration(value=kb_utp.concentration, units='M')
        model_ppi_c.concentration = wc_lang.Concentration(value=kb_ppi.concentration, units='M')
        model_h2o_c.concentration = wc_lang.Concentration(value=kb_h2o.concentration, units='M')
        model_h_c.concentration = wc_lang.Concentration(value=kb_h.concentration, units='M')

        # get or create RNA species        
        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna in rnas:
            species_type = model.species_types.get_or_create(id=rna.id)
            species = species_type.species.get_or_create(compartment=cytosol)
            species.concentration = wc_lang.Concentration(value=rna.concentration, units='M')

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        submodel = self.submodel
        cytosol = model.compartments.get_one(id='c')
        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        cell = self.knowledge_base.cell
        kb_rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for kb_rna in kb_rnas:
            rxn = submodel.reactions.get_or_create(id=kb_rna.id.replace('rna_', 'transcription_'))

            model_rna = model.species_types.get_one(id=kb_rna.id).species.get_one(compartment=cytosol)
            seq = kb_rna.get_seq()
            rxn.participants = []
            rxn.participants.create(species=atp, coefficient=-seq.count('A'))
            rxn.participants.create(species=ctp, coefficient=-seq.count('C'))
            rxn.participants.create(species=gtp, coefficient=-seq.count('G'))
            rxn.participants.create(species=utp, coefficient=-seq.count('U'))
            rxn.participants.create(species=h, coefficient=-(kb_rna.get_len() - 1))
            rxn.participants.create(species=model_rna, coefficient=1)
            rxn.participants.create(species=ppi, coefficient=kb_rna.get_len())
            rxn.participants.create(species=h2o, coefficient=kb_rna.get_len() - 1)

    def gen_rate_laws(self):
        """ Generate rate laws associated with submodel """
        cell = self.knowledge_base.cell

        mean_volume = cell.properties.get_one(id='mean_volume').value
        mean_doubling_time = cell.properties.get_one(id='mean_doubling_time').value

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna, rxn in zip(rnas, self.submodel.reactions):
            rl = rxn.rate_laws.create()
            rl.direction = wc_lang.RateLawDirection.forward
            rl.equation = wc_lang.RateLawEquation(expression='k_cat')
            rl.k_cat = rna.concentration * scipy.constants.Avogadro * mean_volume * numpy.log(2) * (1 / rna.half_life)
            rl.k_m = None
