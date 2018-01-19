""" Make random sequence of polymers

:Author: Cathy Wang <cathy_wang@college.harvard.edu>
:Date: 2018-01-14
:Copyright: 2018, Karr Lab
:License: MIT
"""

import os
import numpy as np

class RandomSeqGen(object):
    """ Make random genetic sequences
    """

    # genetic constants
    DNA_NUCLEOTIDES = ['A','T','C','G']
    DNA_START_CODON = 'TAC'
    DNA_STOP_CODON = 'ATC'
    NUCLEOTIDE_COMPLEMENT = {'A':'U','T':'A','C':'G','G':'C'}
    DEFAULT_MIN_PROT_LEN = 100
    DEFAULT_MAX_PROT_LEN = 1000
    RNA_CODON_TRANS_FILE = os.path.join(os.path.dirname(__file__), 'data/fixtures/codon_translation.csv')

    def __init__(self):

        # read and store codon translation table
        self.translation_table = {}

        for line in open(RandomSeqGen.RNA_CODON_TRANS_FILE):
            # use '#' to identify header or other comment lines, partly disregarding RFC 4180
            if line[0] == '#': continue
            line   = line.rstrip('\n')
            (rna_codon, amino_acid, molecular_weight, charge) = line.split(',')

            # key = each codon in translation table
            # values = corresponding amino acid, molecular weight, and charge
            self.translation_table[rna_codon] = [amino_acid, float(molecular_weight), int(charge)]

    def make_genome(self, prot_lengths):
        """ Make random genome (DNA nucleotide) sequence

        Args:
            prot_lengths (:obj:`list`): list of lengths of each protein

        Returns:
            :obj:`string`: genome sequence for cell
        """
        genome = ''
        for i in range(len(prot_lengths)):
            genome += RandomSeqGen.DNA_START_CODON

            for j in range(prot_lengths[i]-2):
                valid = False
                while not(valid):
                    nuc1 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    nuc2 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    nuc3 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    codon = RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc1] + RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc2] + RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc3]

                    if codon in self.translation_table.keys():
                        valid = True
                genome += nuc1 + nuc2 + nuc3

            genome += RandomSeqGen.DNA_STOP_CODON
        return genome

    def rna_transcript(self, genome):
        """ Complement genome to find mRNA transcript

        Args:
            genome (:obj:`string`): genome (DNA nucleotide) sequence for cell (includes all genes)

        Returns:
            :obj:`string`: mRNA nucleotide sequence for cell
        """
        rna = ''
        for nuc in genome:
            rna += RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc[i]]

        return rna

    def prot_data(self, protein):
        """ Find amino acid sequence, molecular weight, and charge of protein

        Args:
            protein (:obj:`string`): nucleotide sequence of protein

        Returns:
            :obj:`string`: amino acid sequence of protein
            :obj:`float`: molecular weight of protein
            :obj:`int`: charge of protein
        """
        aminoacid = ''
        molweight = 0
        charge = 0

        for i in range(0,len(protein),3):
            codon = protein[i:i+3]

            if codon in self.translation_table.keys():
                aminoacid += self.translation_table[codon][0]
                molweight += self.translation_table[codon][1]
                charge += self.translation_table[codon][2]

        return aminoacid, molweight, charge

    def prot_len(self, num_genes, min_prot_len=DEFAULT_MIN_PROT_LEN, max_prot_len=DEFAULT_MAX_PROT_LEN):
        """ Randomly generate lengths (in amino acids) of all proteins in cell

        Args:
            num_genes (:obj:`int`): number of genes in cell
            min_prot_len (:obj:`int`, optional): shortest desired protein length
            max_prot_len (:obj:`int`, optional): longest desired protein length

        Returns:
            :obj:`list`: list of ints for length of each protein
        """
        prot_length = np.random.randint(min_prot_len, max_prot_len, num_genes)
        return prot_length

    def make_genes(self, genome):
        """ Parse genome to create list of DNA nucleotide (gene) sequences for the proteins in cell

        Args:
            genome (:obj:`string`): genome (DNA nucleotide) sequence for cell (includes all genes)

        Returns:
            :obj:`list`: list of DNA nucleotide sequences of each protein
        """
        genes = ['']
        codon = ''
        index = 0

        for j in range(0, len(genome), 3):
            codon = genome[j:j+3]
            genes[index] += codon
            if codon == RandomSeqGen.DNA_STOP_CODON:
                index += 1
                genes.append('')

        # get rid of last empty string in list
        genes = genes[:-1]

        return genes

    def make_mrnas(self, genes):
        """ Create list of RNA nucleotide sequences for the proteins in cell

        Args:
            genes (:obj:'list'): list of DNA nucleotide sequences of each protein

        Returns:
            :obj:`list`: list of RNA nucleotide sequences of each protein
        """
        prot_rna = []

        for nuc in genes:
            mrna = ''
            for i in range(len(nuc)):
                mrna += RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc[i]]

            prot_rna.append(mrna)

        return prot_rna

    def make_proteins(self, prot_rna):
        """ Make lists of data (amino acid sequence, molecular weight, and charge) for the proteins in cell

        Args:
            prot_rna (:obj:`list`): list of mRNA nucleotide sequences of each protein

        Returns:
            :obj:`tuple`: lists of amino acid sequences, molecular weights, and charges corresponding to each protein
        """
        prot_aa = []
        prot_mw = []
        prot_charge = []

        for i in range(len(prot_rna)):
            aminoacid, molweight, charge = self.prot_data(prot_rna[i])

            prot_aa.append(aminoacid)
            prot_mw.append(molweight)
            prot_charge.append(charge)

        return (prot_aa, prot_mw, prot_charge)

    def gen_species_types(self, num_genes):
        """ Creates and stores genetic and structural data for each protein

        Args:
            num_genes (:obj:`int`): number of genes in cell

        Returns:
            :obj:`tuple`: lists for DNA, RNA, and amino acid sequence/molecular weight/charge data for each protein
        """
        prot_lengths = self.prot_len(num_genes)
        genome = self.make_genome(prot_lengths)
        genes = self.make_genes(genome)
        rnas = self.make_mrnas(genes)
        proteins = self.make_proteins(rnas)

        return (genes, rnas, proteins)

    def output_sequences(self, proteins, out_file):
        """ Make compiled lists of data (amino acid sequence, molecular weight, and charge) for the proteins in cell

        Args:
            proteins (:obj:`tuple`): lists of amino acid sequences, molecular weights, and charges corresponding to each protein
            out_file (:obj:`string`): destination file for writing data

        """
        proteins = np.array(proteins)
        proteins = proteins.transpose()

        with open(out_file, 'w') as outfile:
            outfile.write('aa sequence , molecular weight , charge \n')

            for p in proteins:
                outfile.write('%s' % p[0] + ',' + '%s' % p[1] + ',' + '%s' % p[2] +'\n')

