"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-19
:Copyright: 2018, Karr Lab
:License: MIT
"""
import math
import scipy.stats as stats
import scipy.constants
import wc_kb
import wc_kb_gen
import numpy
import math
from numpy import random
from Bio.Seq import Seq, Alphabet
from Bio.Data import CodonTable

class ObservablesGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates observable objects for proteins and tRNAs that are assigned to specific functions. Adds these observables to the knowledge base.

    Options:
    
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        options = self.options
        assigned_trnas = options.get('assigned_trnas', ['tRNA_Ser', 'tRNA_Leu', 'tRNA_Arg',
                                                        'tRNA_Thr', 'tRNA_Gly', 'tRNA_Phe',
                                                        'tRNA_Trp', 'tRNA_Lys', 'tRNA_Ile',
                                                        'tRNA_Ala', 'tRNA_Met', 'tRNA_Gln',
                                                        'tRNA_Pro', 'tRNA_Val', 'tRNA_Cys',
                                                        'tRNA_Tyr', 'tRNA_His', 'tRNA_Asn',
                                                        'tRNA_Asp', 'tRNA_Glu'])

        rnas = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.RnaSpeciesType)

        trnas = []
        for rna in rnas:
            if rna.type == wc_kb.RnaType.tRna:
                trnas.append(rna)
        
        assert (len(assigned_trnas) <= len(trnas))
        options['assigned_trnas'] = assigned_trnas

        assigned_proteins = options.get('assigned_proteins', ['IF1', 'IF2', 'IF3', 'EFtu', 'EFts',
                                                              'EFg', 'RF1', 'RF2', 'RF3',
                                                              'deg_ATPase', 'deg_protease', 'deg_rnase'])

        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.ProteinSpeciesType)

        assert(len(assigned_proteins) <= len(prots))
        options['assigned_proteins'] = assigned_proteins



    def gen_components(self):
        self.assign_species()
        self.assign_observables()


    def assign_observables(self):
        assigned_trnas = self.options['assigned_trnas']
        assigned_proteins = self.options['assigned_proteins']
        cell = self.knowledge_base.cell

        for trna in assigned_trnas:
            kb_trna = cell.species_types.get(id = trna)
            obs_trna = cell.observables.get_or_create(id = kb_trna.id)
            obs_trna.name = kb_trna.name
            

        for prot in assigned_proteins:
            kb_prot = cell.species_types.get(id = prot)
            obs_prot = cell.observables.get_or_create(id = kb_prot.id)
            obs_prot.name = kb_prot.name


   def assign_species(self):
        """ Takes random samples of the generated rnas and proteins and assigns them functions based on the included list of proteins and rnas"""

        assigned_trnas = self.options['assigned_trnas']
        assigned_proteins = self.options['assigned_proteins']

        prots = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.ProteinSpeciesType)
        rnas = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.RnaSpeciesType)

        trnas = []
        for rna in rnas:
            if rna.type == wc_kb.RnaType.tRna:
                trnas.append(rna)

        sampled_trnas = numpy.random.choice(
            trnas, len(assigned_trnas), replace=False)

        assigned_trnas = iter(assigned_trnas)

        for rna in sampled_trnas:
            rna_name = next(assigned_trnas)
            rna.id = rna_name
            rna.name = rna_name

        sampled_proteins = numpy.random.choice(
            prots, len(assigned_proteins), replace=False)

        assigned_proteins = iter(assigned_proteins)
        for protein in sampled_proteins:
            protein_name = next(assigned_proteins)
            protein.id = protein_name
            protein.name = protein_name
                

   
