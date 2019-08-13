""" Test the generation of a random chromosome and accompanying mRNA and protein sequences
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""
import math
import wc_kb
import os
import shutil
import tempfile
import unittest
from rand_wc_model_gen.kb_gen import genome
from Bio.Data import CodonTable
from rand_wc_model_gen import kb_gen

@unittest.skip("broken_legacy")
class TestGenomeGenerator(unittest.TestCase):

    def setUp(self):

        self.dir = tempfile.mkdtemp()

        # The knowledge base that is needed for the KBComponentGenerator
        # Creates the GenomeGenerator object and sets the parameters as given
        self.whole_gen = kb_gen.KbGenerator(options={
            'component': {
                'PropertiesGenerator': {
                    'mean_volume': 1e-15,
                    'mean_doubling_time': 1000.,
                },
                'GenomeGenerator': {
                    'mean_num_genes': 200.,
                    'seq_path': os.path.join(self.dir, 'kb_seq.fna'),
                },
            },
        })

        self.kb = self.whole_gen.run()

        component_options = self.whole_gen.options.get('component', {})

        self.options = component_options.get('GenomeGenerator', {})

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_init(self):
        self.assertEqual(type(self.whole_gen), kb_gen.KbGenerator)

    def test_num_chromosomes(self):
        chromosomes = self.kb.cell.species_types.get(
            __type=wc_kb.core.DnaSpeciesType)
        self.assertEqual(len(chromosomes),
                         self.options.get('num_chromosomes'))

    def test_rna_props(self):
        rRna = 0
        tRna = 0
        sRna = 0
        rnas = self.kb.cell.species_types.get(
            __type=wc_kb.prokaryote_schema.RnaSpeciesType)

        for rna in rnas:
            if rna.type == wc_kb.core.RnaType.rRna:
                rRna += 1
            elif rna.type == wc_kb.core.RnaType.tRna:
                tRna += 1
            elif rna.type == wc_kb.core.RnaType.sRna:
                sRna += 1

        total = len(rnas)

        rRna_prop = rRna / total
        tRna_prop = tRna / total
        sRna_prop = sRna / total

        real_rRna = self.options.get('rRNA_prop')
        real_tRna = self.options.get('tRNA_prop')
        real_sRna = self.options.get('ncRNA_prop')

        self.assertAlmostEqual(
            rRna_prop, real_rRna, delta=3 * math.sqrt(real_rRna))
        self.assertAlmostEqual(
            tRna_prop, real_tRna, delta=3 * math.sqrt(real_tRna))
        self.assertAlmostEqual(
            sRna_prop, real_sRna, delta=3 * math.sqrt(real_sRna))

    # test total number of RNAs (should match number of transcription units)
    # test total number of proteins (should match number of GeneLocus objects with mRNA)
    def test_rna_num(self):
        rnas = self.kb.cell.species_types.get(
            __type=wc_kb.prokaryote_schema.RnaSpeciesType)
        tus = self.kb.cell.loci.get(
            __type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
        self.assertEqual(len(rnas), len(tus))

    def test_prot_num(self):
        prots = self.kb.cell.species_types.get(
            __type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        genes = self.kb.cell.loci.get(
            __type=wc_kb.prokaryote_schema.GeneLocus)
        geneCount = 0
        for gene in genes:
            if gene.type == wc_kb.core.GeneType.mRna:
                geneCount += 1
        self.assertEqual(geneCount, len(prots))

    def test_start_codon(self):
        trans_table = self.kb.translation_table
        START_CODONS = CodonTable.unambiguous_dna_by_id[trans_table].start_codons

        genes = self.kb.cell.loci.get(__type=wc_kb.prokaryote_schema.GeneLocus)
        for gene in genes:
            if gene.type == wc_kb.core.GeneType.mRna:

                self.assertIn(gene.get_seq()[0:3], START_CODONS)

    def test_stop_codon(self):
        trans_table = self.kb.translation_table
        genes = self.kb.cell.loci.get(__type=wc_kb.prokaryote_schema.GeneLocus)
        STOP_CODONS = CodonTable.unambiguous_dna_by_id[trans_table].stop_codons
        for gene in genes:
            if gene.type == wc_kb.core.GeneType.mRna:
                self.assertIn(gene.get_seq()[-3:], STOP_CODONS)

    # tests first and last sites of proteins
    def test_prot_sites(self):
        pass

    def test_length(self):
        # Tests that the average length of the genes is within 3 standard deviations of the expected.

        genes = self.kb.cell.loci.get(__type=wc_kb.prokaryote_schema.GeneLocus)

        sum_len = 0
        for gene in genes:
            sum_len += gene.get_len()
        avg_len = (sum_len / len(genes)) / 3
        mean_gen_len = self.options.get('mean_gene_len')

        self.assertAlmostEqual(
            avg_len, mean_gen_len, delta=3 * math.sqrt(mean_gen_len))

        # checks average lengths of 5'/3' UTRs on transcription units with mRna
    def test_utrs(self):
        tus = self.kb.cell.loci.get(
            __type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
        sum_five_prime = 0
        sum_three_prime = 0
        mRnaCount = 0
        for tu in tus:
            if tu.genes[0].type == wc_kb.core.GeneType.mRna:  # checks if it is mRna
                mRnaCount += 1
                five_prime_gene = tu.genes[0]
                three_prime_gene = tu.genes[len(tu.genes)-1]
                sum_five_prime += abs(five_prime_gene.start - tu.start)
                sum_three_prime += abs(three_prime_gene.end - tu.end)

        five_prime_len = self.options.get('five_prime_len')
        three_prime_len = self.options.get('three_prime_len')

        self.assertAlmostEqual(sum_five_prime/mRnaCount,
                               five_prime_len, delta=3 * math.sqrt(five_prime_len))
        self.assertAlmostEqual(sum_three_prime/mRnaCount,
                               three_prime_len, delta=3 * math.sqrt(three_prime_len))

    def test_operons(self):
        tus = self.kb.cell.loci.get(
            __type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
        gene_sum = 0
        operonCount = 0
        for tu in tus:
            if len(tu.genes) > 1:  # if operon
                operonCount += 1
                gene_sum += len(tu.genes)

        avg_operon_gen = gene_sum / operonCount
        operon_gen_num = self.options.get('operon_gen_num')

        genes = self.kb.cell.loci.get(__type=wc_kb.prokaryote_schema.GeneLocus)
        geneCount = 0

        for gene in genes:
            if gene.type == wc_kb.core.GeneType.mRna:
                geneCount += 1

        avg_in_operon = gene_sum / geneCount
        operon_prop = self.options.get('operon_prop')

        self.assertAlmostEqual(avg_in_operon, operon_prop,

                               delta=3 * math.sqrt(operon_prop))
        self.assertAlmostEqual(
            avg_operon_gen, operon_gen_num, delta=3 * math.sqrt(operon_gen_num))

    def test_protein_start_codon(self):

        for protein in self.kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType):
            seq = str(protein.get_seq())
            self.assertEqual(seq[0], 'M')

    def test_assignment(self):
        rna = self.kb.cell.species_types.get(name='tRNA_TTT')
        rna = rna[0]
        assert (rna.type == wc_kb.core.RnaType.tRna)

        protein = self.kb.cell.species_types.get(id="translation_init_factors")
        assert(type(protein[0]) == wc_kb.prokaryote_schema.ProteinSpeciesType)
