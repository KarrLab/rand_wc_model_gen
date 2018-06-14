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

GEN_LEN = 10
INTER_LEN = 10
GEN_NUM = 20
TRANSLATION_TABLE = 1


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
        self.indexlist = self.gen.indexList

        # generates a list of the genes and of the intergenic regions
        self.genes = []
        self.intergenes = []

        lastend = None

        for tuple in self.indexlist:
            self.genes.append(self.seq_str[tuple[0]:tuple[1]+1])

            if lastend:
                self.intergenes.append(self.seq_str[lastend:tuple[0]])
            lastend = tuple[1]

    def test_run(self):
        self.assertIsInstance(self.seq, Seq)

    def test_length(self):

        # Tests that the average length of the genes + intergenenic sequences is within 3 standard deviations of the expected.
        # TODO: check if this is an appropriate test, change if not
        self.assertAlmostEqual(
            len(self.seq)/GEN_NUM, (GEN_LEN+INTER_LEN)*3, delta=5*math.sqrt((INTER_LEN/GEN_NUM+GEN_LEN/GEN_NUM)))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
