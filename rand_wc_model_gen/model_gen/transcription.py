""" Generator for models based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import wc_model_gen


class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for transcription submodel """

    def gen_compartments(self):
        model = self.model
        cytosol = model.compartments.get_or_create(id='c')
        cytosol.name = 'cytosol'

    def gen_species(self):
        """ Generate species associated with submodel """
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

        # get or create metabolite species
        atp = model.species_types.get_or_create(id='atp')
        ctp = model.species_types.get_or_create(id='ctp')
        gtp = model.species_types.get_or_create(id='gtp')
        utp = model.species_types.get_or_create(id='utp')
        ppi = model.species_types.get_or_create(id='ppi')
        h2o = model.species_types.get_or_create(id='h2o')
        h = model.species_types.get_or_create(id='h')

        atp.species.get_or_create(compartment=cytosol)
        ctp.species.get_or_create(compartment=cytosol)
        gtp.species.get_or_create(compartment=cytosol)
        utp.species.get_or_create(compartment=cytosol)
        ppi.species.get_or_create(compartment=cytosol)
        h2o.species.get_or_create(compartment=cytosol)
        h.species.get_or_create(compartment=cytosol)

        # get or create RNA species
        for gene in self.knowledge_base.cell.loci.get(__type=wc_kb.GeneLocus):
            species_type = model.species_types.get_or_create(id=gene.id.replace('gene_', 'rna_'))
            species_type.species.get_or_create(compartment=cytosol)

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

        for gene in self.knowledge_base.cell.loci.get(__type=wc_kb.GeneLocus):
            rxn = submodel.reactions.create(id=gene.id.replace('gene_', 'transcription_'))

            rna = model.species_types.get_one(id=gene.id.replace('gene_', 'rna_')).species.get_one(compartment=cytosol)
            seq = gene.get_seq().transcribe()
            rxn.participants.create(species=atp, coefficient=-seq.count('A'))
            rxn.participants.create(species=ctp, coefficient=-seq.count('C'))
            rxn.participants.create(species=gtp, coefficient=-seq.count('G'))
            rxn.participants.create(species=utp, coefficient=-seq.count('U'))
            rxn.participants.create(species=h, coefficient=-(gene.get_len() - 1))
            rxn.participants.create(species=rna, coefficient=1)
            rxn.participants.create(species=ppi, coefficient=gene.get_len())
            rxn.participants.create(species=h2o, coefficient=gene.get_len() - 1)

    def gen_rate_laws(self):
        """ Generate rate laws associated with submodel """
        pass  # pragma: no cover
