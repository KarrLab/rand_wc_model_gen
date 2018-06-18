""" Test the generation of a random chromosome and accompanying mRNA and protein sequences
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb
import unittest
from rand_wc_model_gen.kb_gen import genome
from Bio.Seq import Seq

num_chromosomes = 1
chromosome_topology = 'circular'
mean_gc_frac = .58
mean_num_genes = 4500
mean_gene_len = 900


class TestGenomeGenerator(unittest.TestCase):

    def setUp(self):
        # The knowledge base that is needed for the KBComponentGenerator
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()

        # Creates the GenomeGenerator object and sets the parameters as given
        self.gen = genome.GenomeGenerator(kb)
        options = self.gen.options
        options['mean_gene_len'] = mean_gene_len
        options['mean_num_genes'] = mean_num_genes

        # Generates the sequences
        self.gen.gen_components()

        # the gene sequence that was generated
        self.seq = self.gen.knowledge_base.cell.species_types[0].get_seq()
        # represented as a string
        self.seq_str = str(self.seq)

    def test_run(self):  # confirm that the program runs and generates a sequence
        self.assertIsInstance(self.seq, Seq)

    def test_start_codon(self):
        for gene in self.genes:
            # print(gene)
            self.assertIn(gene[0:3], self.gen.START_CODONS)

    def test_stop_codon(self):
        for gene in self.genes:
            self.assertIn(gene[-3:], self.gen.STOP_CODONS)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
