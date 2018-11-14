"""
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
::Date: 2018-07-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from Bio.Data import CodonTable
import numpy
import wc_kb
import wc_kb_gen


class ObservablesGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates observable objects for proteins and tRNAs that are assigned to specific functions. Adds these observables to the knowledge base.

    Options:
        * assigned_trnas (:obj:`list`): A list of the names of trnas to be created
        * assigned_proteins (:obj:`list`): A list of the names of proteins to be created
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        options = self.options
        default_assigned_trnas = []

        translation_table = int(options.get('translation_table', 1))
        codon_table = CodonTable.unambiguous_dna_by_id[translation_table]
        for codon in codon_table.forward_table.keys():
            default_assigned_trnas.append('tRNA_{}'.format(codon))
        assigned_trnas = options.get('assigned_trnas', default_assigned_trnas)

        rnas = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.prokaryote_schema.RnaSpeciesType)

        count = 0
        for rna in rnas:
            if rna.type == wc_kb.core.RnaType.tRna:
                count += 1

        assert (len(assigned_trnas) <= count)
        options['assigned_trnas'] = assigned_trnas

        assigned_proteins = options.get('assigned_proteins', [
            'rna_polymerase',

            'ribosome',
            'translation_init_factors',
            'translation_elongation_factors',
            'translation_release_factors',
            
            'deg_atpase',
            'degrade_rnase',
            'degrade_protease',
            ])

        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.prokaryote_schema.ProteinSpeciesType)

        assert(len(assigned_proteins) <= len(prots))
        options['assigned_proteins'] = assigned_proteins

    def gen_components(self):
        """ Takes random samples of the generated rnas and proteins and assigns them functions based on the included list of proteins and rnas"""
        cell = self.knowledge_base.cell
        cytosol = cell.compartments.get_one(id='c')

        assigned_trnas = self.options['assigned_trnas']
        assigned_proteins = self.options['assigned_proteins']

        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        rnas = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.prokaryote_schema.RnaSpeciesType)

        trnas = []
        for rna in rnas:
            if rna.type == wc_kb.core.RnaType.tRna:
                trnas.append(rna)

        sampled_trnas = numpy.random.choice(
            trnas, len(assigned_trnas), replace=False)

        assigned_trnas = iter(assigned_trnas)

        for rna in sampled_trnas:
            rna_name = next(assigned_trnas)
            rna.id = rna_name
            rna.name = rna_name
            observable = cell.observables.get_or_create(id=rna_name+'_obs')
            observable.name = rna_name
            #print(observable.name)
            species = rna.species.get_or_create(compartment=cytosol)
            species_coeff = species.species_coefficients.get_or_create(coefficient=1)
            observable.species.append(species_coeff)

        sampled_proteins = numpy.random.choice(
            prots, len(assigned_proteins), replace=False)

        assigned_proteins = iter(assigned_proteins)
        for protein in sampled_proteins:
            protein_name = next(assigned_proteins)
            protein.id = protein_name
            protein.name = protein_name
            observable = cell.observables.get_or_create(id=protein_name+'_obs')
            observable.name = protein_name
            species = protein.species.get_or_create(compartment=cytosol)
            species_coeff = species.species_coefficients.get_or_create(coefficient=1)
            observable.species.append(species_coeff)
