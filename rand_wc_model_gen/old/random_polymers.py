""" Make random sequence of polymers

:Author: Cathy Wang <cathy_wang@college.harvard.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-01-14
:Copyright: 2018, Karr Lab
:License: MIT
"""

import os
import numpy as np

# TODO(Arthur): double-check that Start and Stop codons are being handled correctly
# TODO(Arthur): create config file
# use it to replace some constants like DEFAULT_MIN_PROT_LEN, DEFAULT_MAX_PROT_LEN, RNA_CODON_TRANS_FILE
# TODO(Arthur): combine this with enrich_polymers.py in core.py
# TODO(Arthur): enable config of PRNG seed for reproducibility
# TODO(Arthur): catch exceptions, such as in write_output()


class RandomSeqGen(object):
    """ Make random genetic sequences

    Attributes:
        translation_table (:obj:`dict`): map of amino acid encoding codons, from codon to
            corresponding AA, molecular weight, and charge
    """

    # genetic constants
    DNA_NUCLEOTIDES = ['A','T','C','G']
    DNA_START_CODON = 'TAC'
    DNA_STOP_CODON = 'ATC'
    NUCLEOTIDE_COMPLEMENT = {'A':'U','T':'A','C':'G','G':'C'}
    PROT_DATA_HEADER = 'aa sequence,molecular weight,charge\n'
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

        The genome consists of `len(prot_lengths)` protein coding genes. Each gene is a sequence of
        valid coding codons, bounded by start and stop codons. All valid coding codons are equally
        likely.

        Args:
            prot_lengths (:obj:`list` of :obj:`int`): the nucleotide lengths of proteins in the genome

        Returns:
            :obj:`string`: DNA nucleotides in genome sequence for cell
        """
        genome = ''
        for i in range(len(prot_lengths)):
            genome += RandomSeqGen.DNA_START_CODON

            for j in range(prot_lengths[i]):
                valid = False
                while not(valid):
                    nuc1 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    nuc2 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    nuc3 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    codon = RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc1] + \
                        RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc2] + RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc3]

                    if codon in self.translation_table.keys():
                        valid = True
                genome += nuc1 + nuc2 + nuc3

            genome += RandomSeqGen.DNA_STOP_CODON
        return genome

    def rna_transcript(self, genome):
        """ Complement genome to find mRNA transcript

        Args:
            genome (:obj:`string`): genome (DNA nucleotide) sequence for cell

        Returns:
            :obj:`string`: mRNA nucleotide sequence for cell
        """
        rna = ''
        for nuc in genome:
            rna += RandomSeqGen.NUCLEOTIDE_COMPLEMENT[nuc]

        return rna

    def prot_data(self, protein):
        """ Find amino acid sequence, molecular weight, and charge of a protein

        Args:
            protein (:obj:`string`): DNA nucleotide sequence of a protein

        Returns:
            (`tuple`):

                (:obj:`string`: amino acid sequence of protein,
                :obj:`float`: molecular weight of protein,
                :obj:`int`: charge of protein)
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
        """ Generate random lengths (in amino acids) for all proteins in a cell

        Args:
            num_genes (:obj:`int`): number of genes in cell
            min_prot_len (:obj:`int`, optional): the shortest desired protein length
            max_prot_len (:obj:`int`, optional): the longest desired protein length

        Returns:
            :obj:`list`: list of ints for length of each protein
        """
        prot_length = np.random.randint(min_prot_len, max_prot_len, num_genes)
        return prot_length

    def make_genes(self, genome):
        """ Parse genome to create list of DNA nucleotide (gene) sequences for the proteins in cell

        Args:
            genome (:obj:`string`): genome (DNA nucleotide) sequence for cell

        Returns:
            :obj:`list`: DNA nucleotide sequence of each protein in `genome`
        """
        genes = []
        gene = ''

        for j in range(0, len(genome), 3):
            codon = genome[j:j+3]
            gene += codon
            if codon == RandomSeqGen.DNA_STOP_CODON:
                genes.append(gene)
                gene = []

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

    def make_proteins(self, mrna_seqs):
        """ Make lists of the properties of the proteins in a cell

        The properties are amino acid sequence, molecular weight, and charge.

        Args:
            mrna_seqs (:obj:`list` of :obj:`str`): mRNA nucleotide sequence for each protein

        Returns:
            (`tuple`): parallel lists of protein properties:

                (:obj:`list` of :obj:`str`: amino acid sequences, :obj:`list` of :obj:`float`: molecular weights,
                :obj:`list` of :obj:`float`: charges)
        """
        prot_aa = []
        prot_mw = []
        prot_charge = []

        for i in range(len(mrna_seqs)):
            aminoacid, molweight, charge = self.prot_data(mrna_seqs[i])

            prot_aa.append(aminoacid)
            prot_mw.append(molweight)
            prot_charge.append(charge)

        return (prot_aa, prot_mw, prot_charge)

    def gen_polymers(self, num_genes):
        """ Generates genetic sequences and structural data for a cell

        Args:
            num_genes (:obj:`int`): number of genes in the cell

        Returns:
            (:obj:`tuple`): the genome, genes, mRNAs, and proteins in a cell; see respective `make_<polymer>()` methods for details

                (:obj:`str`: DNA nucleotides, :obj:`tuple` of :obj:`list`: gene properties,
                :obj:`tuple` of :obj:`list`: mRNA properties,
                :obj:`tuple` of :obj:`list`: protein properties)
        """
        prot_lengths = self.prot_len(num_genes)
        genome = self.make_genome(prot_lengths)
        genes = self.make_genes(genome)
        mrnas = self.make_mrnas(genes)
        proteins = self.make_proteins(mrnas)

        return (genome, genes, mrnas, proteins)

    def format_prot_data(self, proteins):
        """ Format data (amino acid sequence, molecular weight, and charge) for the proteins in cell

        Args:
            proteins (:obj:`tuple`): lists of amino acid sequences, molecular weights, and charges
                corresponding to each protein

        Return:
            prot_data_string (:obj:`list`): strings of structural data corresponding to each protein

        """
        proteins = np.array(proteins)
        proteins = proteins.transpose()
        prot_data_string = []

        for p in proteins:
            prot_data_string.append('%s' % p[0] + ',' + '%s' % p[1] + ',' + '%s' % p[2] +'\n')

        return prot_data_string

    def write_output(self, prot_data_string, outfile, header=PROT_DATA_HEADER):
        """ Write protein data to a file

        Args:
            prot_data_string (:obj:`list`): strings of amino acid sequences, molecular weights, and
                charges corresponding to each protein
            out_file (:obj:`string`): destination file for writing data
            header (:obj:`string`, optional): titles for column of data in `out_file`
        """
        with open(outfile, 'w') as outfile:
            outfile.write(header)
            for s in prot_data_string:
                outfile.write(s)
