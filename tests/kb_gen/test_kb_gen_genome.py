""" Test the generation of a random chromosome and accompanying mRNA and protein sequences
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""
import math
import wc_kb
import unittest
from rand_wc_model_gen.kb_gen import genome


class TestGenomeGenerator(unittest.TestCase):

    def setUp(self):
        # The knowledge base that is needed for the KBComponentGenerator
        kb = wc_kb.KnowledgeBase()

        # Creates the GenomeGenerator object and sets the parameters as given
        options = {
            'mean_num_genes': 1000}
        self.gen = genome.GenomeGenerator(kb, options)
        self.gen.run()

    def test_init(self):
        self.assertEqual(type(self.gen), genome.GenomeGenerator)

    def test_num_chromosomes(self):
        chromosomes = self.gen.knowledge_base.cell.species_types.get(
            __type=wc_kb.core.DnaSpeciesType)
        self.assertEqual(len(chromosomes),
                         self.gen.options.get('num_chromosomes'))

    def test_rna_props(self):
        rRna = 0
        tRna = 0
        sRna = 0
        rnas = self.gen.knowledge_base.cell.species_types.get(
            __type=wc_kb.core.RnaSpeciesType)

        for rna in rnas:
            if rna.type == wc_kb.RnaType.rRna:
                rRna += 1
            elif rna.type == wc_kb.RnaType.tRna:
                tRna += 1
            elif rna.type == wc_kb.RnaType.sRna:
                sRna += 1

        total = len(rnas)

        rRna_prop = rRna / total
        tRna_prop = tRna / total
        sRna_prop = sRna / total

        real_rRna = self.gen.options.get('rRNA_prop')
        real_tRna = self.gen.options.get('tRNA_prop')
        real_sRna = self.gen.options.get('ncRNA_prop')

        self.assertAlmostEqual(
            rRna_prop, real_rRna, delta=3 * math.sqrt(real_rRna))
        self.assertAlmostEqual(
            tRna_prop, real_tRna, delta=3 * math.sqrt(real_tRna))
        self.assertAlmostEqual(
            sRna_prop, real_sRna, delta=3 * math.sqrt(real_sRna))

    # test total number of RNAs (should match number of transcription units)
    # test total number of proteins (should match number of GeneLocus objects with mRNA)
    def test_rna_num(self):
        rnas = self.gen.knowledge_base.cell.species_types.get(
            __type=wc_kb.core.RnaSpeciesType)
        tus = self.gen.knowledge_base.cell.loci.get(
            __type=wc_kb.core.TranscriptionUnitLocus)
        self.assertEqual(len(rnas), len(tus))

    def test_prot_num(self):
        prots = self.gen.knowledge_base.cell.species_types.get(
            __type=wc_kb.core.ProteinSpeciesType)
        genes = self.gen.knowledge_base.cell.loci.get(
            __type=wc_kb.core.GeneLocus)
        geneCount = 0
        for gene in genes:
            if gene.type == wc_kb.core.GeneType.mRna:
                geneCount += 1
        self.assertEqual(geneCount, len(prots))

    def test_start_codon(self):
        trans_table = self.gen.knowledge_base.translation_table
        START_CODONS = trans_table.start_codons
        STOP_CODONS = trans_table.stop_codons
        genes = self.gen.knowledge_base.cell.loci.get(__type=wc_kb.GeneLocus)
        for gene in genes:
            if gene.type == wc_kb.GeneType.mRna:

                self.assertIn(gene.get_seq()[0:3], START_CODONS)

    def test_stop_codon(self):
        trans_table = self.gen.knowledge_base.translation_table
        genes = self.gen.knowledge_base.cell.loci.get(__type=wc_kb.GeneLocus)
        START_CODONS = trans_table.start_codons
        STOP_CODONS = trans_table.stop_codons
        for gene in genes:
            if gene.type == wc_kb.GeneType.mRna:
                self.assertIn(gene.get_seq()[-3:], STOP_CODONS)

    # tests first and last sites of proteins
    def test_prot_sites(self):
        pass

    def test_length(self):
        # Tests that the average length of the genes is within 3 standard deviations of the expected.

        genes = self.gen.knowledge_base.cell.loci.get(__type=wc_kb.GeneLocus)

        sum_len = 0
        for gene in genes:
            sum_len += gene.get_len()
        avg_len = (sum_len / len(genes)) / 3
        mean_gen_len = self.gen.options.get('mean_gene_len')

        self.assertAlmostEqual(
            avg_len, mean_gen_len, delta=3 * math.sqrt(mean_gen_len))

        # checks average lengths of 5'/3' UTRs on transcription units with mRna
    def test_utrs(self):
        tus = self.gen.knowledge_base.cell.loci.get(
            __type=wc_kb.core.TranscriptionUnitLocus)
        sum_five_prime = 0
        sum_three_prime = 0
        mRnaCount = 0
        for tu in tus:
            if tu.genes[0].type == wc_kb.GeneType.mRna:  # checks if it is mRna
                mRnaCount += 1
                five_prime_gene = tu.genes[0]
                three_prime_gene = tu.genes[len(tu.genes)-1]
                sum_five_prime += abs(five_prime_gene.start - tu.start)
                sum_three_prime += abs(three_prime_gene.end - tu.end)

        five_prime_len = self.gen.options.get('five_prime_len')
        three_prime_len = self.gen.options.get('three_prime_len')

        self.assertAlmostEqual(sum_five_prime/mRnaCount,
                               five_prime_len, delta=3 * math.sqrt(five_prime_len))
        self.assertAlmostEqual(sum_three_prime/mRnaCount,
                               three_prime_len, delta=3 * math.sqrt(three_prime_len))

    def test_operons(self):
        tus = self.gen.knowledge_base.cell.loci.get(
            __type=wc_kb.core.TranscriptionUnitLocus)
        gene_sum = 0
        operonCount = 0
        for tu in tus:
            if len(tu.genes) > 1:  # if operon
                operonCount += 1
                gene_sum += len(tu.genes)

        avg_operon_gen = gene_sum / operonCount
        operon_gen_num = self.gen.options.get('operon_gen_num')

        genes = self.gen.knowledge_base.cell.loci.get(__type=wc_kb.GeneLocus)
        geneCount = 0

        for gene in genes:
            if gene.type == wc_kb.GeneType.mRna:
                geneCount += 1

        avg_in_operon = gene_sum / geneCount
        operon_prop = self.gen.options.get('operon_prop')

        self.assertAlmostEqual(avg_in_operon, operon_prop,

                               delta=3 * math.sqrt(operon_prop))
        self.assertAlmostEqual(
            avg_operon_gen, operon_gen_num, delta=3 * math.sqrt(operon_gen_num))

    '''def test_protein_start_codon(self):
        for protein in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            seq = str(protein.get_seq())
            self.assertEqual(seq[0], 'M')

    def test_protein_stop_codon(self):
        for protein in self.gen.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            seq = str(protein.get_seq())
            self.assertEqual(seq[-1], '*')'''

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
