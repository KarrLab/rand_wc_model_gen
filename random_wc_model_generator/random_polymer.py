""" Make random sequence of polymers

:Author: Cathy Wang <cathy_wang@college.harvard.edu>
:Date: 2018-01-14
:Copyright: 2018, Karr Lab
:License: MIT
"""

import os
import numpy as np

class RandomSeqGen(object):
    """ Make random sequences
    """

    # genetic constants
    RNA_NUCLEOTIDES = ['A','U','C','G']
    START_CODON = 'AUG'
    STOP_CODON = 'UAG'
    NUCLEOTIDE_COMP = {'A':'T','U':'A','C':'G','G':'C'}
    CODON_TRANS_FILE = os.path.join(os.path.dirname(__file__), 'data/fixtures/codon_translation.csv')

    def __init__(self):
    
        # read and store codon translation table
        self.translation_table = {}

        for line in open(CODON_TRANS_FILE):
            if line[0] == '#': continue      
            line   = line.rstrip('\n')       

            fields = line.split(',')
    
            # key = each codon in translation table
            # values = corresponding amino acid, molecular weight, and charge
            self.translation_table[fields[0]] = fields[1:4]

    def prot_seq(self, length):
        """ Randomly create nucleotide sequence for protein 

        Args:
            length (:obj:`int`): desired length (in amino acids) of protein

        Returns:
            :obj:`string`: nucleotide sequence of protein
        """
        protein = random_seq_gen.START_CODON

        for i in range(length-2):
            valid = False
        
            while not(valid):
                nuc1 = '{}'.format(np.random.choice(random_seq_gen.RNA_NUCLEOTIDES))
                nuc2 = '{}'.format(np.random.choice(random_seq_gen.RNA_NUCLEOTIDES))
                nuc3 = '{}'.format(np.random.choice(random_seq_gen.RNA_NUCLEOTIDES))
                codon = nuc1 + nuc2 + nuc3
                
                if codon in self.translation_table.keys():
                    valid = True
            protein += codon
    
        protein += random_seq_gen.STOP_CODON

        return protein
    
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
                molweight += float(self.translation_table[codon][1])
                charge += int(self.translation_table[codon][2])

        return aminoacid, molweight, charge
    
    def dna_seq(self, rna_sequence):
        """ Find DNA sequence

        Args:
            rna_sequence (:obj:`string`): RNA nucleotide sequence

        Returns:
            :obj:`string`: DNA sequence of cell
        """
        dna_sequence = ''
    
        for nuc in rna_sequence:
            dna_sequence += random_seq_gen.NUCLEOTIDE_COMP[nuc]
    
        return dna_sequence

    def prot_len(self, num_genes):
        """ Randomly generate lengths (in amino acids) of all proteins in cell

        Args:
            num_genes (:obj:`int`): number of genes in cell

        Returns:
            :obj:`list`: list of ints for length of each protein
        """
        prot_len = np.random.randint(50,1000,num_genes)
        return prot_len
    
    def make_genes(self, num_genes, prot_len):
        """ Create compiled list of nucleotide sequences for the proteins in cell

        Args:
            num_genes (:obj:`int`): number of genes in cell
            prot_len (:obj:'list'): list of ints for length of each protein

        Returns:
            :obj:`list`: list of nucleotide sequences of each protein
        """
        genes = []
        for i in range(num_genes):
            prot_nucl_seq = self.prot_seq(prot_len[i])
            genes.append(prot_nucl_seq)
        return genes
    
    def make_rna(self, genes):
        """ Find RNA sequence

        Args:
            genes (:obj:`list`): list of nucleotide sequences of all proteins

        Returns:
            :obj:`string`: RNA sequence of cell
        """
        rna_sequence = ''.join(genes)

        return rna_sequence

    def make_proteins(self, genes):
        """ Make compiled lists of data (amino acid sequence, molecular weight, and charge) for the proteins in cell

        Args:
            genes (:obj:`list`): list of nucleotide sequences of each protein

        Returns:
            :obj:`tuple`: lists of amino acid sequences, molecular weights, and charges corresponding to each protein
        """
        prot_aa = []
        prot_mw = []
        prot_charge = []
        
        for i in range(len(genes)):
            aminoacid, molweight, charge = self.prot_data(genes[i])
            
            prot_aa.append(aminoacid)
            prot_mw.append(molweight)
            prot_charge.append(charge)
            
        return (prot_aa, prot_mw, prot_charge)

    def gen_species_types(self, num_genes):
        """ Creates and stores genetic and structural data for each protein

        Args:
            num_genes (:obj:`int`): number of genes in cell

        Returns:
            :obj:`tuple`: lists for nucleotide sequences of each protein and RNA and amino acid sequence/molecular weight/charge data
        """
        prot_lengths = self.prot_len(num_genes)
        genes = self.make_genes(num_genes,prot_lengths)
        rna = self.make_rna(genes)
        proteins = self.make_proteins(genes)
        return (genes, rna, proteins)

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


random_seq_gen = RandomSeqGen()
NUM_GENES = 10

(genes, rna, proteins) = random_seq_gen.gen_species_types(NUM_GENES)
OUT_FILE = 'protein_data.csv'
random_seq_gen.output_sequences(proteins, OUT_FILE)
