""" Test polymer sequence generation

:Author: Cathy Wang <cathy_wang@college.harvard.edu>
:Date: 2018-01-17
:Copyright: 2018, Karr Lab
:License: MIT
"""
import unittest
import numpy as np

from random_wc_model_generator.random_polymer import RandomSeqGen

NUM_GENES=10
NUM_RUNS=100

class TestGenerateSeq(unittest.TestCase):
    
    def setUp(self):
        self.rsg = RandomSeqGen()

    def test_length(self):
        ''' Compare length of protein (in amino acids) with length assigned by randomized protein length list
        '''
        prot_lengths = self.rsg.prot_len(num_genes=NUM_GENES)
        genome = self.rsg.make_genome(prot_lengths)
        genes = self.rsg.make_genes(genome)
        rnas = self.rsg.make_mrnas(genes)
        proteins = self.rsg.make_proteins(rnas)
        
        for i in range(len(rnas)):
            # account for lack of stop codon amino acid with +1
            self.assertEqual(len(proteins[0][i])+1, prot_lengths[i])
        
    def test_base_stats(self):
        ''' Check abundance of each base in gene nucleotide sequence
        '''        
        base_abundance = np.zeros((NUM_RUNS,4))
        for i in range(NUM_RUNS):
            prot_length = self.rsg.prot_len(num_genes=NUM_GENES)
            genome = self.rsg.make_genome(prot_length)

            base_abundance[i][0] = genome.count('A')/len(genome)
            base_abundance[i][1] = genome.count('C')/len(genome)
            base_abundance[i][2] = genome.count('T')/len(genome)
            base_abundance[i][3] = genome.count('G')/len(genome)
        
        for j in range(0,4):
            np.testing.assert_allclose(np.mean(base_abundance[:][j]), 0.25, rtol = 0.2)

    def test_prot_data_reading(self):
        ''' Check valid values are being read and stored from protein data file
        '''
        (genes, rnas, proteins) = self.rsg.gen_species_types(NUM_GENES)
        for i in range(NUM_GENES):
            aa, mw, charge = self.rsg.prot_data(rnas[i])
            self.assertGreater(mw,0)
            self.assertIsInstance(charge, int)
