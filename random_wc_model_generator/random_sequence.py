""" desc

:Author: 
:Date: 
:Copyright: 2018, Karr Lab
:License: MIT
"""

import os
import np

'''
constants in UPPER_CASE
organize code in classes
more descriptive names
docstrings
PEP8/Google coding style
    lower_case, UPPER_CASE, CamelCase names
testing code

'''

class RandomSeqGen(object):
    """ Make random sequences
    """

    # genetic constants
    RNA_NUCLEOTIDES = ['A','U','C','G']
    START_CODON = 1
    STOP_CODON = 1
    NUCLEOTIDE_COMP = 2
    CODON_TRANS_FILE = os.path.join(os.path.dirname(__file__), 'data/fixtures/codon_translation.csv')
    
    # read and store codon translation table
    # codon_trans_file = 'translationTable.csv'

    def __init__(self):
    
        self.translation_table = {}

        for line in open(CODON_TRANS_FILE):
            if line[0] == '#': continue      
            line   = line.rstrip('\n')       

            fields = line.split(',')
    
            # key = each codon in translation table, value = corresponding amino acid
            self.translation_table[fields[0]] = fields[1]

    def prot_seq(self, length):
        """ Parse the command line

        Args:
            length (:obj:`int`): desired length (in amino acids) of protein

        Returns:
            :obj:`string`: nucleotide sequence for each protein
        """
        # start codon
        protein = 'AUG'

        for i in range(length-2):
            valid = False
        
            while not(valid):
                nuc1 = '{}'.format(np.random.choice(random_seq_gen.RNA_NUCLEOTIDES))
                nuc2 = '{}'.format(np.random.choice(random_seq_gen.RNA_NUCLEOTIDES))
                nuc3 = '{}'.format(np.random.choice(random_seq_gen.RNA_NUCLEOTIDES))
                codon = nuc1 + nuc2 + nuc3
                # if codon exits in list of valid codons (translationTable)
                if codon in self.translation_table.keys():
                    valid = True
            protein += codon
    
        # stop codon
        protein += 'UAG'

        return protein

    def RNAseq(self, protein_sequence):

        rna_sequence = ''.join(protein_sequence)

        return rna_sequence

    # complement of RNA sequence
    def DNAseq(self, rna_sequence):
        dna_sequence = ''
        complement = {'A':'T','U':'A','C':'G','G':'C'}
    
        for nuc in rna_sequence:
            dna_sequence += complement[nuc]
    
        return dna_sequence

    protein_sequence = []
    # number of proteins
    n = 3
    # create list of random lengths (in terms of amino acids) for each of the n proteins
    # could be potentially imported, with set length of genes
    protLen = np.random.randint(50,1000,n)

    for i in range(n):
        protein = protSeq(protLen[i])
        # list of nucleotide sequences for each protein
        protein_sequence.append(protein)

    rna_sequence = RNAseq(protein_sequence)
    dna_sequence = DNAseq(rna_sequence)

    def make_genes(self, num_genes):
        pass

    def gen_species_types(self, num_genes, output_file):
        genes = self.make_genes(num_genes)
        rnas = self.make_rna(genes)
        proteins = self.make_protein(rnas)
        return (genes, rnas, proteins)

    def output_sequences(genes, rnas, proteins, OUT_FILE):
        pass

random_seq_gen = RandomSeqGen()
NUM_GENES = 10
(genes, rnas, proteins) = random_seq_gen.gen_species_types(NUM_GENES)
OUT_FILE = 'test.tsv'
random_seq_gen.output_sequences(genes, rnas, proteins, OUT_FILE)
