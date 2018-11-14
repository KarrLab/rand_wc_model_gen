""" Tests of generation of chromosomes and genes

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen.kb_gen import chrs_genes_tus
import numpy
import unittest
import wc_kb


class ChromosomesGenesTusGeneratorTestCase(unittest.TestCase):
    def test_rand(self):
        gen = chrs_genes_tus.ChromosomesGenesTusGenerator(None)
        ints = gen.rand(1000, 10000)
        numpy.testing.assert_equal(numpy.ceil(ints), ints)
        self.assertGreater((numpy.mean(ints)-1000) / 1000, -1e-2)
        self.assertLess((numpy.mean(ints)-1000) / 1000, 1e-2)

    def test_run(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        mean_num_genes = 200        
        gen = chrs_genes_tus.ChromosomesGenesTusGenerator(kb, options={
            'num_chromosomes': 2,
            'mean_gc_frac': 1,
            'mean_num_genes': mean_num_genes,
            'mean_gene_len': 100,
            'mean_coding_frac': 0.75,
        })
        gen.run()

        self.assertEqual(len(cell.species_types), 2)

        gc = 0
        chr_len = 0
        for chr in cell.species_types:
            gc += chr.seq.count('C') + chr.seq.count('G')
            chr_len += len(chr.seq)
        self.assertEqual(gc / chr_len, 1)

        tus = cell.loci.get(__type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
        genes = cell.loci.get(__type=wc_kb.prokaryote_schema.GeneLocus)
        self.assertAlmostEqual(len(tus), mean_num_genes, delta=5 * numpy.sqrt(mean_num_genes))
        self.assertAlmostEqual(len(genes), mean_num_genes, delta=5 * numpy.sqrt(mean_num_genes))

        mu = 100 / 0.75 * 0.25 / 2
        self.assertAlmostEqual(genes[0].start, mu, delta=5 * numpy.sqrt(mu))

        gene_len = 0
        gene_pos = 0
        for gene in genes:
            gene_len += gene.get_len()
            gene_pos += gene.strand == wc_kb.core.PolymerStrand.positive
        self.assertAlmostEqual(gene_len / len(genes), 100, delta=5 * numpy.sqrt(100/100))
        self.assertAlmostEqual(gene_pos / len(genes), 0.5, delta=5 * numpy.sqrt((1 - 0.5) * 0.5 / 100))

        self.assertAlmostEqual(gene_len / chr_len, 0.75, delta=5 * numpy.sqrt((1 - 0.75) * 0.75 / 100))
