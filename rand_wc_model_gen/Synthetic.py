"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-06
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import math
import numpy as np
import random
from Bio.Seq import Seq

class Synthetic(object):
    """
    Creates synthetic chromosome with randomized genes/intergenic regions. Creates RNA and protein objects corresponding to the genes on this chromosome. Associates the chromosome, RNAs, proteins
    with a knowledge base object (and its Cell attribute)
    
    """

    def __init__(self, gen_len, inter_len, gen_num):
        """
        Args:
            gen_len (:obj:`int`): average gene length
            inter_len (:obj:`int`): average intergenic region length
            gen_num (:obj:`int`): number of genes on chromosome

        """
        self.kb = wc_kb.KnowledgeBase() #knowledge base instance has attribute cell
        cell = wc_kb.Cell()
        self.kb.cell = cell
        self.kb.translation_table = 1 #translation table for RNA codons

        indexList = self.create_chromosome(gen_len, inter_len, gen_num) #indexList of start/end positions of each gene, creates 'synthetic' chromosome

        self.create_rnas_proteins(gen_num, indexList) #creates RNA and protein objects corresponding to the genes on chromosome


    def create_chromosome(self, gen_len, inter_len, gen_num):
        """ Creates 'synthetic' chromsome with randomized genes/intergenic regions

        Args:
            gen_len (:obj:`int`): average gene length
            inter_len (:obj:`int`): average intergenic region length
            gen_num (:obj:`int`): number of genes on chromosome


        Returns:
            :obj:`list`: list of tuples of start and end positions of each gene on chromosome
        """
        gene_dist = np.random.normal(gen_len, math.sqrt(gen_len/gen_num), gen_num).tolist() #takes random samples out of Gaussian distribution with mean of average gene length
        gene_dist = [round(x) for x in gene_dist]
        inter_dist = np.random.normal(inter_len, math.sqrt(inter_len/gen_num), gen_num).tolist() #takes random samples out of Gaussian distribution with mean of average intergenic length
        inter_dist = [round(x) for x in inter_dist]
        chromosome = wc_kb.DnaSpeciesType()
        seq = ''
        arr = ['A', 'G', 'C', 'T']
        START_CODON = 'ATG' #start codon
        STOP_CODONS = ['TAG', 'TAA', 'TGA'] #stop codons
        indexList = []
        index = 1
        for i in range(2 * gen_num):
            if i % 2 == 0: #if i is even, region is a gene
                gene = ''
                gen_length = gene_dist.pop()
                indexList.append((index, index + (gen_length * 3) - 1))
                gene += START_CODON #add start codon
                for k in range(gen_length - 2):
                    codon = 'TAG'
                    while codon in STOP_CODONS: #to make sure random codon is not stop codon (premature termination)
                       codon = random.choice(arr) + random.choice(arr) + random.choice(arr) #create random codon
                    gene += codon #add randomly chosen codon
                gene += random.choice(STOP_CODONS) #add stop codon
                index += gen_length * 3
                seq += gene
                
            else: #if i is odd, region is intergenic
                inter_length = inter_dist.pop()
                index += inter_length * 3
                for k in range(inter_length):
                    seq += random.choice(arr) + random.choice(arr) + random.choice(arr) #add randomly chosen base triple

        
        chromosome.seq = Seq(seq) #associate the random chromosome sequence with the DnaSpeciesType object
        self.kb.cell.species_types.append(chromosome) #add chromosome to kb.cell speciestypes list

        return indexList


    def create_rnas_proteins(self, gen_num, indexList):
        """ Creates RNA and protein objects corresponding to genes on chromosome

        Args:
            gen_num (:obj:`int`): number of genes on chromosome
            indexList (:obj: 'list'): list of tuples of start and end positions of each gene on chromosome

        """

        chromosome = self.kb.cell.species_types[0]

        for i in range(gen_num):
            rna = wc_kb.RnaSpeciesType() #creates RnaSpeciesType for RNA sequence corresponding to gene
            gene = wc_kb.GeneLocus() #GeneLocus object for gene sequence, attribute of ProteinSpeciesType object
            gene.start = indexList[i][0] 
            gene.end = indexList[i][1]
            gene.polymer = chromosome
            gene.cell = self.kb.cell
            transcription_locus = wc_kb.TranscriptionUnitLocus() #TranscriptionUnitLocus object - attribute of RnaSpeciesType object, associated with gene sequence
            transcription_locus.polymer = chromosome
            transcription_locus.start = gene.start
            transcription_locus.end = gene.end

            rna.transcription_units.append(transcription_locus) 
            
            self.kb.cell.species_types.append(rna) #adds corresponding mRNA sequence to speciestypes list of kb.cell

            prot = wc_kb.ProteinSpeciesType() #creates ProteinSpeciesType object for corresponding protein sequence

            prot.gene = gene #associates protein with GeneLocus object for corresponding gene
            
            self.kb.cell.species_types.append(prot) #adds ProteinSpeciesType object to kb.cell speciestypes list
                



        
