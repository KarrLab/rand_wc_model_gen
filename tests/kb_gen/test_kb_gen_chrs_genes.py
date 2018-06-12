""" Tests of generation of chromosomes and genes

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen.kb_gen import chrs_genes
import numpy
import unittest
import wc_kb


class ChromosomesGenesGeneratorTestCase(unittest.TestCase):
    def test_rand(self):
        gen = chrs_genes.ChromosomesGenesGenerator(None)
        ints = gen.rand(1000, 10000)
        numpy.testing.assert_equal(numpy.ceil(ints), ints)
        self.assertGreater((numpy.mean(ints)-1000) / 1000, -1e-2)
        self.assertLess((numpy.mean(ints)-1000) / 1000, 1e-2)

    def test_gen_components(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        gen = chrs_genes.ChromosomesGenesGenerator(kb, options={
            'num_chromosomes': 2,
            'avg_gc_frac': 1,
            'avg_num_genes': 100,
            'avg_gene_len': 100,
            'avg_coding_frac': 0.75,
        })
        gen.gen_components()

        self.assertEqual(len(cell.species_types), 2)

        gc = 0
        chr_len = 0
        for chr in cell.species_types:
            gc += chr.seq.count('C') + chr.seq.count('G')
            chr_len += len(chr.seq)
        self.assertEqual(gc / chr_len, 1)

        self.assertAlmostEqual(len(cell.loci), 100, delta=25)

        self.assertAlmostEqual(cell.loci[0].start, 100 / 0.75 * 0.25 / 2, delta=20)

        gene_len = 0
        for gene in cell.loci:
            gene_len += gene.get_len()
        self.assertAlmostEqual(gene_len / len(cell.loci), 100, delta=5)

        self.assertAlmostEqual(gene_len / chr_len, 0.75, delta=0.1)
