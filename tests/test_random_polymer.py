""" Test polymer sequence generation

:Author: Cathy Wang <cathy_wang@college.harvard.edu>
:Date: 2018-01-17
:Copyright: 2018, Karr Lab
:License: MIT
"""
import unittest

from random_polymer import RandomSeqGen as rsg

class TestGenerateSeq(unittest.TestCase):
    NUM_GENES=10
    
    def test_length(self):
        ''' Compare length of protein (in amino acids) with length assigned by randomized protein length list
        '''
        prot_len = rsg.prot_len(num_genes=NUM_GENES)
        genes = rsg.make_genes(num_genes=NUM_GENES, prot_len)
        proteins = rsg.make_proteins(genes)
        
        for i in range(len(genes)):
            self.assertEqual(len(proteins[0][i]), prot_len[i])
        
    def test_base_stats(self):
        ''' Check abundance of each base in gene nucleotide sequence
        '''
        
        prot_length = rsg.prot_len(num_genes=NUM_GENES)
        genes = rsg.make_genes(num_genes=NUM_GENES, prot_len)
        
        base_abundance = np.zeros(len(genes),4)
        for i in range(len(genes)):
            base_abundance[i][0] = genes[i].count('A')/len(genes[i])
            base_abundance[i][1] = genes[i].coun('C')/len(genes[i])
            base_abundance[i][2] = genes[i].count('U')/len(genes[i])
            base_abundance[i][3] = genes[i].count('G')/len(genes[i])
        for j in range(0,4):
            np.testing.assert_allclose(mean(base_abundance[:][j]), 0.25, rtol = 0.2)

    def test_prot_data_reading(self):
        ''' Check valid values are being read and stored from protein data file
        '''
        protein = rsg.prot_seq(length=100)
        aa, mw, charge = rsg.prot_data(protein)
        self.assertGreater(mw,0)
        self.assertIsInstance(charge, int)

        
