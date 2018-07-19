"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-06-06
:Copyright: 2018, Karr Lab
:License: MIT
"""
import math
import scipy.stats as stats
import scipy.constants
import wc_kb
import wc_kb_gen
import numpy
import math
from numpy import random
from Bio.Seq import Seq, Alphabet
from Bio.Data import CodonTable


class GenomeGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates synthetic chromosome with randomized genes/intergenic regions. Creates RNA and protein objects corresponding to the genes on this chromosome. Associates the chromosome, RNAs, proteins
    with a knowledge base object (and its Cell attribute)

    Options:

    * num_chromosomes (:obj:`int`): number of chromosomes
    * mean_gc_frac (:obj:`float`): fraction of nucleotides which are G or C
    * mean_num_genes (:obj:`float`): mean number of genes
    * mean_gene_len (:obj:`float`): mean codon length of a gene
    * mean_coding_frac (:obj:`float`): mean coding fraction of the genome
    * translation_table (:obj:'int'): The NCBI standard genetic code used
    * ncRNA_prop (:obj:'float'): The proportion of non coding RNAs
    * rRNA_prop  (:obj:'float'): The proportion of ribosomal RNAs
    * tRNA_prop (:obj:'float'): The proportion of transfer RNAs
    * five_prime_len (:obj:'int'): Average 5' UTR length for transcription units
    * three_prime_len (:obj:'int'): Average 3' UTR length for transcription units
    * operon_prop (:obj:'float'): Proportion of genes that should be in an operon (polycistronic mRNA)
    * operon_gen_num (:obj:'int'): Average number of genes in an operon
    * mean_copy_number (:obj:`float`): mean copy number of each RNA
    * mean_half_life (:obj:`float`): mean half-life of RNAs
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        # Default options are loosely  based on Escherichia coli K-12
        # Nucleic Acids Research 41:D605-12 2013

        options = self.options

        num_chromosomes = options.get('num_chromosomes', 1)
        assert(num_chromosomes >= 1 and int(
            num_chromosomes) == num_chromosomes)
        options['num_chromosomes'] = num_chromosomes

        chromosome_topology = options.get('chromosome_topology', 'circular')
        assert(chromosome_topology in ['circular', 'linear'])
        options['chromosome_topology'] = chromosome_topology

        mean_gc_frac = options.get('mean_gc_frac', 0.58)
        assert(mean_gc_frac >= 0 and mean_gc_frac <= 1)
        options['mean_gc_frac'] = mean_gc_frac

        mean_num_genes = options.get('mean_num_genes', 4500)
        assert(mean_num_genes >= 1)
        options['mean_num_genes'] = mean_num_genes

        # print(mean_num_genes)

        ncRNA_prop = options.get('ncRNA_prop', 0.014)
        assert(ncRNA_prop >= 0 and ncRNA_prop <= 1)
        options['ncRNA_prop'] = ncRNA_prop

        rRNA_prop = options.get('rRNA_prop', 0.0056)
        assert(rRNA_prop >= 0 and rRNA_prop <= 1)
        options['rRNA_prop'] = rRNA_prop

        tRNA_prop = options.get('tRNA_prop', 0.5)
        assert(tRNA_prop >= 0 and tRNA_prop <= 1)
        options['tRNA_prop'] = tRNA_prop

        assert((ncRNA_prop + rRNA_prop + tRNA_prop) <= 1)

        # DOI: 10.1093/molbev/msk019
        mean_gene_len = options.get(
            'mean_gene_len', 308)  # codon length (924 bp)
        assert(mean_gene_len >= 1)
        options['mean_gene_len'] = mean_gene_len

        # DOI: 10.1007/s10142-015-0433-4
        mean_coding_frac = options.get('mean_coding_frac', 0.88)
        assert(mean_coding_frac > 0 and mean_coding_frac < 1)
        options['mean_coding_frac'] = mean_coding_frac

        translation_table = int(options.get('translation_table', 1))
        assert(translation_table in range(1, 32))
        options['translation_table'] = translation_table

        five_prime_len = int(options.get('five_prime_len', 7))
        assert(five_prime_len >= 0)
        options['five_prime_len'] = five_prime_len

        three_prime_len = int(options.get('three_prime_len', 5))  # guess
        assert(three_prime_len >= 0)
        options['three_prime_len'] = three_prime_len

        operon_prop = (options.get('operon_prop', 0.2))  # guess
        assert(operon_prop >= 0 and operon_prop <= 1)
        options['operon_prop'] = operon_prop

        operon_gen_num = int(options.get('operon_gen_num', 3))
        assert(operon_gen_num >= 2)
        options['operon_gen_num'] = operon_gen_num

        # DOI: 10.1038/ismej.2012.94
        mean_copy_number = options.get('mean_copy_number', 0.4)
        assert(mean_copy_number > 0)
        options['mean_copy_number'] = mean_copy_number

        # DOI: 10.1073/pnas.0308747101
        mean_half_life = options.get('mean_half_life', 2.1 * 60)
        assert(mean_half_life > 0)
        options['mean_half_life'] = mean_half_life

    def gen_components(self):
        self.gen_genome()
        self.gen_tus()
        self.gen_rnas_proteins()

    def gen_genome(self):
        '''Construct knowledge base components and generate the DNA sequence'''

        # get options
        options = self.options
        num_chromosomes = options.get('num_chromosomes')
        mean_gene_len = options.get('mean_gene_len')
        translation_table = options.get('translation_table')
        mean_num_genes = options.get('mean_num_genes')
        mean_coding_frac = options.get('mean_coding_frac')
        mean_gc_frac = options.get('mean_gc_frac')
        chromosome_topology = options.get('chromosome_topology')
        ncRNA_prop = options.get('ncRNA_prop')
        rRNA_prop = options.get('rRNA_prop')
        tRNA_prop = options.get('tRNA_prop')

        # print(tRNA_prop)

        cell = self.knowledge_base.cell

        self.knowledge_base.translation_table = translation_table

        codon_table = CodonTable.unambiguous_dna_by_id[
            translation_table]

        # start codons from NCBI list
        START_CODONS = codon_table.start_codons

        # stop codons from NCBI list
        STOP_CODONS = codon_table.stop_codons

        BASES = ['A', 'C', 'G', 'T']  # The DNA Bases

        # The probability of each base being selected randomly
        PROB_BASES = [(1 - mean_gc_frac) / 2, mean_gc_frac /
                      2, mean_gc_frac/2, (1-mean_gc_frac)/2]

        # Create a chromosome n times
        for i_chr in range(num_chromosomes):
            # number of genes in the chromosome
            num_genes = self.rand(mean_num_genes / num_chromosomes)[0]
            # list of gene lengths (generated randomly) on chromosome
            gene_lens = 3 * self.rand(mean_gene_len, count=num_genes, min=2)

            intergene_lens = 3 * self.rand(
                mean_gene_len / mean_coding_frac * (1 - mean_coding_frac), count=num_genes)

            # sequence base triple length
            seq_len = numpy.sum(gene_lens) + numpy.sum(intergene_lens)

            seq_str = []
            # generates seq based on random codons (NOT start/stop codons)
            for i in range(0, seq_len, 3):
                codon_i = STOP_CODONS[0]

                codon_i = "".join(random.choice(
                    BASES, p=PROB_BASES, size=(3,)))

                seq_str.append(codon_i)

            seq_str = "".join(seq_str)
            seq = Seq(seq_str, Alphabet.DNAAlphabet())

            chro = cell.species_types.get_or_create(
                id='chr_{}'.format(i_chr + 1), __type=wc_kb.DnaSpeciesType)
            chro.name = 'Chromosome {}'.format(i_chr + 1)
            chro.circular = chromosome_topology == 'circular'
            chro.double_stranded = True
            chro.seq = seq

            gene_starts = numpy.int64(numpy.cumsum(numpy.concatenate(([0], gene_lens[0:-1])) +
                                                   numpy.concatenate((numpy.round(intergene_lens[0:1] / 2), intergene_lens[1:]))))

            # creates GeneLocus objects for the genes and labels their GeneType (which type of RNA they transcribe)
            for i_gene, gene_start in enumerate(gene_starts):
                gene = self.knowledge_base.cell.loci.get_or_create(
                    id='gene_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.GeneLocus)
                gene.start = gene_start + 1  # 1-indexed
                gene.polymer = chro
                gene.end = gene.start + gene_lens[i_gene] - 1  # 1-indexed
                # print(gene_lens[i_gene] % 3 == 0)
                gene.name = 'gene {} {}'.format(i_chr+1, i_gene+1)
                typeList = [wc_kb.GeneType.mRna, wc_kb.GeneType.rRna,
                            wc_kb.GeneType.sRna, wc_kb.GeneType.tRna]
                prob_rna = [1 - ncRNA_prop - tRNA_prop -
                            rRNA_prop, rRNA_prop, ncRNA_prop, tRNA_prop]
                gene.type = random.choice(typeList, p=prob_rna)
                if gene.type == wc_kb.GeneType.mRna:  # if mRNA, then set up start/stop codons in the gene
                    start_codon = random.choice(START_CODONS)
                    stop_codon = random.choice(STOP_CODONS)
                    seq_str = str(chro.seq)
                    seq_str = seq_str[:gene.start-1] + start_codon + \
                        seq_str[gene.start+2: gene.end-3] + \
                        stop_codon + seq_str[gene.end:]
                    for i in range(gene.start+2, gene.end-3, 3):
                        # print(seq_str[i:i+3])
                        while seq_str[i:i+3] in START_CODONS or seq_str[i:i+3] in STOP_CODONS:
                            # print('here')
                            codon_i = "".join(random.choice(
                                BASES, p=PROB_BASES, size=(3,)))
                            seq_str = seq_str[:i]+codon_i+seq_str[i+3:]

                    chro.seq = Seq(seq_str, Alphabet.DNAAlphabet())

    def gen_rnas_proteins(self):
        """ Creates RNA and protein objects corresponding to genes on chromosome

        """
        options = self.options
        mean_copy_number = options.get('mean_copy_number')
        mean_half_life = options.get('mean_half_life')
        mean_volume = self.knowledge_base.cell.properties.get_one(
            id='mean_volume').value
        for chromosome in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.DnaSpeciesType):
            for i in range(len(chromosome.loci)):

                locus = chromosome.loci[i]

                if type(locus) == wc_kb.TranscriptionUnitLocus:
                    tu = locus

                    # creates RnaSpeciesType for RNA sequence corresponding to gene
                    rna = self.knowledge_base.cell.species_types.get_or_create(
                        id='rna_{}'.format(tu.id), __type=wc_kb.RnaSpeciesType)
                    rna.name = 'rna {}'.format(tu.id)
                    # GeneLocus object for gene sequence, attribute of ProteinSpeciesType object
                    if tu.genes[0].type == wc_kb.GeneType.mRna:
                        rna.type = wc_kb.RnaType.mRna
                    elif tu.genes[0].type == wc_kb.GeneType.rRna:
                        rna.type = wc_kb.RnaType.rRna
                    elif tu.genes[0].type == wc_kb.GeneType.tRna:
                        rna.type = wc_kb.RnaType.tRna
                    elif tu.genes[0].type == wc_kb.GeneType.sRna:
                        rna.type = wc_kb.RnaType.sRna

                    # print(rna.type)
                    rna.concentration = random.gamma(
                        1, mean_copy_number) / scipy.constants.Avogadro / mean_volume
                    rna.half_life = random.normal(
                        mean_half_life, numpy.sqrt(mean_half_life))

                    rna.transcription_units.append(tu)

                    # print(rna.get_seq

                    if rna.type == wc_kb.RnaType.mRna:
                        for gene in tu.genes:
                            # creates ProteinSpecipe object for corresponding protein sequence(s)
                            # print(gene.get_seq()[0:3])
                            prot = self.knowledge_base.cell.species_types.get_or_create(
                                id='prot_{}'.format(gene.id), __type=wc_kb.ProteinSpeciesType)
                            prot.name = 'prot_{}'.format(gene.id)

                            prot.cell = self.knowledge_base.cell
                            prot.cell.knowledge_base = self.knowledge_base

                            prot.gene = gene  # associates protein with GeneLocus object for corresponding gene
                            prot.rna = rna
                            prot.half_life = 1
                            prot.concentration = rna.concentration

    def gen_tus(self):
        """ Creates transcription units with 5'/3' UTRs, polycistronic mRNAs, and other types of RNA (tRNA, rRNA, sRNA)

        """

        # validate these options in the options method
        options = self.options
        # 7 bp default (E. coli, wikipedia)
        five_prime_len = options.get('five_prime_len')
        three_prime_len = options.get('three_prime_len')  # 5 bp default guess
        operon_prop = options.get('operon_prop')  # 0.2 default guess
        # 3 genes default (https://academic.oup.com/gbe/article/5/11/2242/653613)
        operon_gen_num = options.get('operon_gen_num')

        for i_chr, chromosome in enumerate(self.knowledge_base.cell.species_types.get(__type=wc_kb.core.DnaSpeciesType)):
            seq = chromosome.seq
            i_gene = 0
            transcription_loci = []
            # print(chromosome.loci)
            # print(len(chromosome.loci))
            # Todo make this into a proper for loop that deals with repeats/additional loci
            while i_gene < len(chromosome.loci):

                gene = chromosome.loci[i_gene]

                if gene.type == wc_kb.GeneType.mRna:
                    # polycistronic mRNA (multiple GeneLocus objects per TranscriptionUnitLocus)
                   # print("mrna")

                    five_prime = self.rand(five_prime_len)[0]
                    three_prime = self.rand(three_prime_len)[0]
                    operon_prob = random.random()

                    # make an operon (polycistronic mRNA, put multiple genes in one TransUnitLocus)
                    if operon_prob <= operon_prop:
                       # print("ploycistronic")
                        operon_genes = self.rand(operon_gen_num)[0]
                        while operon_genes <= 1:
                            operon_genes = self.rand(operon_gen_num)[0]

                        # add 3', 5' UTRs to the ends of the transcription unit (upstream of first gene, downstream of last gene)
                        tu = self.knowledge_base.cell.loci.get_or_create(
                            id='tu_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.TranscriptionUnitLocus)
                        tu.name = 'tu {} {}'.format(i_chr+1, i_gene+1)

                        five_prime_start = gene.start - five_prime
                        if five_prime_start < 0:
                            five_prime_start = 0
                        tu.genes.append(gene)
                        tu.start = five_prime_start

                       # print(chromosome.loci)
                        # print(operon_genes)
                        for k in range(operon_genes-1):

                            i_gene += 1

                            if i_gene >= len(chromosome.loci):
                                break

                            if (chromosome.loci[i_gene]).type == wc_kb.GeneType.mRna:
                                gene = chromosome.loci[i_gene]
                                tu.genes.append(gene)

                            else:
                                break

                        three_prime_end = gene.end + three_prime
                        if three_prime_end >= len(seq):
                            three_prime_end = len(seq) - 1
                        tu.end = three_prime_end
                        transcription_loci.append(tu)

                    else:  # make an individual transcription unit for the gene
                        # print("single")
                        five_prime_start = gene.start - five_prime
                        three_prime_end = gene.end + three_prime
                        if five_prime_start < 0:
                            five_prime_start = 0
                        if three_prime_end >= len(seq):
                            three_prime_end = len(seq) - 1
                        tu = self.knowledge_base.cell.loci.get_or_create(
                            id='tu_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.TranscriptionUnitLocus)
                        tu.start = five_prime_start
                        tu.end = three_prime_end
                        tu.name = 'tu {} {}'.format(i_chr+1, i_gene+1)
                        tu.genes.append(gene)
                        transcription_loci.append(tu)

                # make a transcription unit that transcribes other types of RNA (tRNA, rRNA, sRNA)
                else:
                   # print("not mrna")
                    tu = self.knowledge_base.cell.loci.get_or_create(
                        id='tu_{}_{}'.format(i_chr + 1, i_gene + 1), __type=wc_kb.TranscriptionUnitLocus)
                    tu.name = 'tu {} {}'.format(i_chr+1, i_gene+1)
                    tu.start = gene.start
                    tu.end = gene.end
                    tu.genes.append(gene)
                    transcription_loci.append(tu)

                i_gene += 1
            for locus in transcription_loci:
                locus.polymer = chromosome

    def rand(self, mean, count=1, min=0, max=numpy.inf):
        """ Generated 1 or more random normally distributed integer(s) with standard deviation equal
        to the square root of the mean value.

        Args:
            mean (:obj:`float`): mean value
            count (:obj:`int`): number of random numbers to generate

        Returns:
            :obj:`int` or :obj:`numpy.ndarray` of :obj:`int`: random normally distributed integer(s)
        """
        a = (min-mean)/numpy.sqrt(mean)
        b = (max - mean)/numpy.sqrt(mean)
        mu = mean

        return numpy.int64(numpy.round(stats.truncnorm.rvs(a, b, loc=mean, scale=numpy.sqrt(mean), size=count)))
