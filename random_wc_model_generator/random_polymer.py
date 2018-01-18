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
    # in DNA
    START_CODON = 'TAC'
    STOP_CODON = 'ATC'
    NUCLEOTIDE_COMP = {'A':'U','T':'A','C':'G','G':'C'}
    MIN_PROT_LEN = 100
    MAX_PROT_LEN = 1000
    CODON_TRANS_FILE = os.path.join(os.path.dirname(__file__), 'data/fixtures/codon_translation.csv')

    def __init__(self):
    
        # read and store codon translation table
        self.translation_table = {}

        for line in open(RandomSeqGen.CODON_TRANS_FILE):
            # use '#' to identify header or other comment lines, partly disregarding RFC 4180
            if line[0] == '#': continue
            line   = line.rstrip('\n')       
            (codon, amino_acid, molecular_weight, charge) = line.split(',')
    
            # key = each codon in translation table
            # values = corresponding amino acid, molecular weight, and charge
            self.translation_table[codon] = [amino_acid, float(molecular_weight), int(charge)]
            
    def make_genome(self, prot_length):
        """ Make random genome (DNA nucleotide) sequence

        Args:
            prot_length (:obj:`list`): list of lengths of each protein

        Returns:
            :obj:`string`: genome sequence for cell
        """
        genome = ''
        for i in range(len(prot_length)):
            genome += RandomSeqGen.START_CODON
            
            for j in range(prot_length[i]-2):
                valid = False
                while not(valid):
                    nuc1 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    nuc2 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    nuc3 = '{}'.format(np.random.choice(RandomSeqGen.DNA_NUCLEOTIDES))
                    codon = RandomSeqGen.NUCLEOTIDE_COMP[nuc1] + RandomSeqGen.NUCLEOTIDE_COMP[nuc2] + RandomSeqGen.NUCLEOTIDE_COMP[nuc3]

                    if codon in self.translation_table.keys():
                        valid = True
                genome += nuc1 + nuc2 + nuc3
                
            genome += RandomSeqGen.STOP_CODON
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
            rna += RandomSeqGen.NUCLEOTIDE_COMP[nuc[i]]
        
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
    
    def prot_len(self, num_genes, min_prot_len=MIN_PROT_LEN, max_prot_len=MAX_PROT_LEN):
        """ Randomly generate lengths (in amino acids) of all proteins in cell

        Args:
            num_genes (:obj:`int`): number of genes in cell
            min_prot_len (:obj:`int`): shortest desired protein length
            min_prot_len (:obj:`int`): longest desired protein length

        Returns:
            :obj:`list`: list of ints for length of each protein
        """
        prot_length = np.random.randint(min_prot_len, max_prot_len, num_genes)
        return prot_length
    
    def make_genes(self, genome):
        """ Parse genome to create compiled list of DNA nucleotide (gene) sequences for the proteins in cell

        Args:
            genome (:obj:`string`): genome (DNA nucleotide) sequence for cell (includes all genes)

        Returns:
            :obj:`list`: list of DNA nucleotide sequences of each protein
        """
        genes = ['']
        codon = ''
        index = 0

        for j in range(index,len(genome),3):
            codon = genome[j:j+3]
            if codon == 'ATC':
                genes[index] += codon
                index += 1
                genes.append('')

            else:
                genes[index] += codon
        
        # get rid of last emtpy string in list
        genes = genes[:-1]
    
        return genes

    def make_mrnas(self, genes):
        """ Create compiled list of RNA nucleotide sequences for the proteins in cell

        Args:
            genes (:obj:'list'): list of DNA nucleotide sequences of each protein 

        Returns:
            :obj:`list`: list of RNA nucleotide sequences of each protein
        """
        prot_rna = []
    
        for nuc in genes:
            mrna = ''
            for i in range(len(nuc)):
                mrna += RandomSeqGen.NUCLEOTIDE_COMP[nuc[i]]
            
            prot_rna.append(mrna)
    
        return prot_rna
    
    def make_proteins(self, prot_rna):
        """ Make compiled lists of data (amino acid sequence, molecular weight, and charge) for the proteins in cell

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

'''
RandomSeqGen = RandomSeqGen()
NUM_GENES = 10
(genes, rnas, proteins) = RandomSeqGen.gen_species_types(NUM_GENES)
OUT_FILE = 'protein_data.csv'
RandomSeqGen.output_sequences(proteins, OUT_FILE)
'''
