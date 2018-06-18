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

GEN_LEN = 20  # the average number of codons in each gene
INTER_LEN = 20  # the average number of codons between genes
GEN_NUM = 500  # the exact number of genes present
TRANSLATION_TABLE = 1  # the codon table to use


class TestGenomeGenerator(unittest.TestCase):

    def setUp(self):
        # The knowledge base that is needed for the KBComponentGenerator
        kb = wc_kb.KnowledgeBase()
        self.cell = cell = kb.cell = wc_kb.Cell()

        # Creates the GenomeGenerator object and sets the parameters as given
        self.gen = genome.GenomeGenerator(kb)
        self.gen.options = {
            'mean_num_genes': 1000
        }

        # Generates the sequences and chromosomes
        self.gen.gen_components()

        self.gen.make_tus()
        
        self.gen.gen_rnas_proteins()

        # generates a list of the genes and of the intergenic regions
        #genes = cell.loci.get(__type=wc_kb.GeneLocus)
        #self.gene_seqs = []
        

        #self.assertEqual(sum([len(seq) for seq in self.gene_seqs]) + sum([len(seq) for seq in self.intergene_seqs]), len(seq_str))


    def test_num_chromosomes(self):
        chromosomes = self.cell.loci.get(__type=wc_kb.DnaSpeciesType)
        self.assertEqual(len(chromosomes), self.gen.options.get('num_chromosomes'))

    def test_num_genes(self):
        genes = 

    def test_rna_props(self):
        rRna = 0
        tRna = 0
        sRna = 0
        rnas = self.cell.loci.get(__type=wc_kb.RnaSpeciesType)
        for rna in rnas:
            elif rna.type == wc_kb.RnaType.rRna:
                rRna += 1
            elif rna.type == wc_kb.RnaType.tRna:
                tRna += 1
            else:
                sRna += 1

        total = len(rnas)

        rRna_prop = rRna / total
        tRna_prop = tRna / total
        sRna_prop = sRna / total


        self.assertAlmostEqual(
            rRna_prop, self.gen.options.get('rRna_prop'), delta=3 * math.sqrt(rRna_prop))
        self.assertAlmostEqual(
            tRna_prop, self.gen.options.get('tRna_prop'), delta=3 * math.sqrt(tRna_prop))
        self.assertAlmostEqual(
            sRna_prop, self.gen.options.get('sRna_prop'), delta=3 * math.sqrt(sRna_prop))

        

        
        
    '''def test_start_codon(self):
        genes = self.cell.loci.get(__type=wc_kb.GeneLocus)
        for gene in genes:
            self.assertIn(gene.get_seq()[0:3], self.gen.START_CODONS)

    def test_stop_codon(self):
        genes = self.cell.loci.get(__type=wc_kb.GeneLocus)
        for gene in genes:
            self.assertIn(gene.get_seq()[-3:], self.gen.STOP_CODONS)

    def test_length(self):
        # Tests that the average length of the genes + intergenenic sequences is within 3 standard deviations of the expected.

        self.assertAlmostEqual(
            len(self.seq) / GEN_NUM, 3 * (GEN_LEN + INTER_LEN), delta=3 * math.sqrt(3 * (INTER_LEN + GEN_LEN) / GEN_NUM))

        genes = self.cell.loci.get(__type=wc_kb.GeneLocus)

        sum_len = 0
        for gene in genes:
            sum_len += gene.get_len()
        avg_len = sum_len / GEN_NUM

        self.assertAlmostEqual(
            avg_len, 3 * GEN_LEN, delta=3 * math.sqrt(3 * GEN_LEN / GEN_NUM))

        sum_len = 0
        for intergene_seq in self.intergene_seqs:
            sum_len += len(intergene_seq)
        avg_len = sum_len / GEN_NUM

        self.assertAlmostEqual(
            avg_len, 3 * INTER_LEN, delta=3 * math.sqrt(3 * INTER_LEN / GEN_NUM))

    def test_rna_start_codon(self):
        for rna in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType):
            seq = str(rna.get_seq())
            self.assertEqual(seq[0:3], 'AUG')

    def test_rna_stop_codon(self):
        STOP_CODONS = ['UAG', 'UAA', 'UGA']
        for rna in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType):
            seq = str(rna.get_seq())
            self.assertIn(seq[-3:], STOP_CODONS)

    def test_protein_start_codon(self):
        for protein in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            seq = str(protein.get_seq())
            self.assertEqual(seq[0], 'M')

    def test_protein_stop_codon(self):
        for protein in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            seq = str(protein.get_seq())
            self.assertEqual(seq[-1], '*')

    def test_number_of_rnas(self):
        count = 0
        for rna in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType):
            count += 1
        
        self.assertEqual(count, GEN_NUM)

    def test_number_of_proteins(self):
        count = 0
        for protein in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            count += 1
        self.assertEqual(count, GEN_NUM)'''

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
