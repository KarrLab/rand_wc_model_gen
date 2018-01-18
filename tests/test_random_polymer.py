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


class TestGenerateSeq(unittest.TestCase):
    
    def setUp(self):
        self.rsg = RandomSeqGen()

    def test_length(self):
        ''' Compare length of protein (in amino acids) with length assigned by randomized protein length list
        '''
        prot_length = self.rsg.prot_len(num_genes=NUM_GENES)
        genes = self.rsg.make_genes(NUM_GENES, prot_length)
        proteins = self.rsg.make_proteins(genes)

        for i in range(len(genes)):
            self.assertEqual(len(proteins[0][i]), prot_length[i])
        
    def test_base_stats(self):
        ''' Check abundance of each base in gene nucleotide sequence
        '''
        
        prot_length = self.rsg.prot_len(num_genes=NUM_GENES)
        genes = self.rsg.make_genes(NUM_GENES, prot_length)
        
        base_abundance = np.zeros((len(genes),4))
        for i in range(len(genes)):
            base_abundance[i][0] = genes[i].count('A')/len(genes[i])
            base_abundance[i][1] = genes[i].count('C')/len(genes[i])
            base_abundance[i][2] = genes[i].count('U')/len(genes[i])
            base_abundance[i][3] = genes[i].count('G')/len(genes[i])
        for j in range(0,4):
            np.testing.assert_allclose(np.mean(base_abundance[:][j]), 0.25, rtol = 0.2)

    def test_prot_data_reading(self):
        ''' Check valid values are being read and stored from protein data file
        '''
        protein = self.rsg.prot_seq(length=100)
        aa, mw, charge = self.rsg.prot_data(protein)
        self.assertGreater(mw,0)
        self.assertIsInstance(charge, int)

        
