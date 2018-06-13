"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-06-06
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import math
import numpy as np
import random
from Bio.Seq import Seq


class GenomeGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates synthetic chromosome with randomized genes/intergenic regions. Creates RNA and protein objects corresponding to the genes on this chromosome. Associates the chromosome, RNAs, proteins
    with a knowledge base object (and its Cell attribute)

    """

    def generate_genome(self, gen_len, inter_len, gen_num, translation_table=1):
        """
        Args:
            gen_len (:obj:`int`): average gene length
            inter_len (:obj:`int`): average intergenic region length
            gen_num (:obj:`int`): number of genes on chromosome
            translation_table (:obj:int'): integer specifying the translation table to use
        """
        self.kb = wc_kb.KnowledgeBase()  # knowledge base instance has attribute cell
        cell = wc_kb.Cell()
        self.kb.cell = cell
        self.kb.translation_table = translation_table  # translation table for RNA codons
        # indexList of start/end positions of each gene, creates 'synthetic' chromosome
        indexList = self.create_chromosome(gen_len, inter_len, gen_num)
        # creates RNA and protein objects corresponding to the genes on chromosome
        self.create_rnas_proteins(gen_num, indexList)

    def create_chromosome(self, gen_len, inter_len, gen_num):
        """ Creates 'synthetic' chromsome with randomized genes/intergenic regions

        Args:
            gen_len (:obj:`int`): average gene length
            inter_len (:obj:`int`): average intergenic region length
            gen_num (:obj:`int`): number of genes on chromosome


        Returns:
            :obj:`list`: list of tuples of start and end positions of each gene on chromosome
        """
        gene_dist = np.random.normal(gen_len, math.sqrt(gen_len/gen_num), gen_num).tolist(
        )  # takes random samples out of Gaussian distribution with mean of average gene length
        gene_dist = [round(x) for x in gene_dist]
        # takes random samples out of Gaussian distribution with mean of average intergenic length
        inter_dist = np.random.normal(inter_len, math.sqrt(
            inter_len/gen_num), gen_num).tolist()
        inter_dist = [round(x) for x in inter_dist]
        chromosome = wc_kb.DnaSpeciesType()
        seq = ''
        arr = ['A', 'G', 'C', 'T']
        START_CODON = 'ATG'  # start codon
        STOP_CODONS = ['TAG', 'TAA', 'TGA']  # stop codons
        indexList = []
        index = 1
        for i in range(2 * gen_num):
            if i % 2 == 0:  # if i is even, region is a gene
                gene = ''
                gen_length = gene_dist.pop()
                indexList.append((index, index + (gen_length * 3) - 1))
                gene += START_CODON  # add start codon
                for k in range(gen_length - 2):
                    codon = 'TAG'
                    # to make sure random codon is not stop codon (premature termination)
                    while codon in STOP_CODONS:
                        # create random codon
                        codon = random.choice(
                            arr) + random.choice(arr) + random.choice(arr)
                    gene += codon  # add randomly chosen codon
                gene += random.choice(STOP_CODONS)  # add stop codon
                index += gen_length * 3
                seq += gene

            else:  # if i is odd, region is intergenic
                inter_length = inter_dist.pop()
                index += inter_length * 3
                for k in range(inter_length):
                    # add randomly chosen base triple
                    seq += random.choice(arr) + \
                        random.choice(arr) + random.choice(arr)

        # associate the random chromosome sequence with the DnaSpeciesType object
        chromosome.seq = Seq(seq)
        # add chromosome to kb.cell speciestypes list
        self.kb.cell.species_types.append(chromosome)

        return indexList

    def create_rnas_proteins(self, gen_num, indexList):
        """ Creates RNA and protein objects corresponding to genes on chromosome

        Args:
            gen_num (:obj:`int`): number of genes on chromosome
            indexList (:obj: 'list'): list of tuples of start and end positions of each gene on chromosome

        """

        chromosome = self.kb.cell.species_types[0]

        for i in range(gen_num):
            # creates RnaSpeciesType for RNA sequence corresponding to gene
            rna = wc_kb.RnaSpeciesType()
            # GeneLocus object for gene sequence, attribute of ProteinSpeciesType object
            gene = wc_kb.GeneLocus()
            gene.start = indexList[i][0]
            gene.end = indexList[i][1]
            gene.polymer = chromosome
            gene.cell = self.kb.cell
            # TranscriptionUnitLocus object - attribute of RnaSpeciesType object, associated with gene sequence
            transcription_locus = wc_kb.TranscriptionUnitLocus()
            transcription_locus.polymer = chromosome
            transcription_locus.start = gene.start
            transcription_locus.end = gene.end

            rna.transcription_units.append(transcription_locus)

            # adds corresponding mRNA sequence to speciestypes list of kb.cell
            self.kb.cell.species_types.append(rna)

            # creates ProteinSpeciesType object for corresponding protein sequence
            prot = wc_kb.ProteinSpeciesType()

            prot.gene = gene  # associates protein with GeneLocus object for corresponding gene

            # adds ProteinSpeciesType object to kb.cell speciestypes list
            self.kb.cell.species_types.append(prot)
