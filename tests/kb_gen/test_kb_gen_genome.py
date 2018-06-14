""" Test the generation of a random chromosome and accompanying mRNA and protein sequences
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""
import math
import wc_kb
import wc_kb_gen
import unittest
import numpy as np
from rand_wc_model_gen.kb_gen import genome
from Bio.Seq import Seq

GEN_LEN = 200  # the average number of codons in each gene
INTER_LEN = 200  # the average number of codons between genes
GEN_NUM = 500  # the exact number of genes present
TRANSLATION_TABLE = 1  # the codon table to use


class TestGenomeGenerator(unittest.TestCase):

    def setUp(self):
        # The knowledge base that is needed for the KBComponentGenerator
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()

        # Creates the GenomeGenerator object and sets the parameters as given
        self.gen = genome.GenomeGenerator(kb)
        self.gen.options = {'gen_len': GEN_LEN, 'inter_len': INTER_LEN,
                            'gen_num': GEN_NUM, 'translation_table': TRANSLATION_TABLE}

        # Generates the sequences
        self.gen.gen_components()

        # the gene sequence that was generated
        self.seq = self.gen.knowledge_base.cell.species_types[0].get_seq()
        # represented as a string
        self.seq_str = str(self.seq)

        # a list of tuples containing the starting and ending indices of each gene
        self.indexlist = self.gen.indexList  # Index starts at 1

        # generates a list of the genes and of the intergenic regions
        self.genes = []
        self.intergenes = []
        lastend = None

        for tuple in self.indexlist:
            self.genes.append(self.seq_str[tuple[0]-1:tuple[1]])

            if lastend:
                self.intergenes.append(self.seq_str[lastend:tuple[0]])
            lastend = tuple[1]-1

       # print(self.genes)
        # print(self.intergenes)

    def test_run(self):
        self.assertIsInstance(self.seq, Seq)

    def test_number_of_genes(self):
        self.assertEqual(len(self.genes), GEN_NUM)
        self.assertEqual(len(self.intergenes), GEN_NUM - 1)

    def test_start_codon(self):
        for gene in self.genes:
            # print(gene)
            self.assertIn(gene[0:3], self.gen.START_CODONS)

    def test_stop_codon(self):
        for gene in self.genes:
            self.assertIn(gene[-3:], self.gen.STOP_CODONS)

    def test_length(self):

        # Tests that the average length of the genes + intergenenic sequences is within 3 standard deviations of the expected.
        # TODO: check if this is an appropriate test, change if not
        self.assertAlmostEqual(
            len(self.seq)/GEN_NUM, (GEN_LEN+INTER_LEN)*3, delta=5*math.sqrt((INTER_LEN/GEN_NUM+GEN_LEN/GEN_NUM)))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
