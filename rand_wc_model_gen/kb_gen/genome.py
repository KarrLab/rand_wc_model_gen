"""
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Ashwin Srinivasan
:Date: 2018-06-25
:Copyright: 2018, Karr Lab
:License: MIT
"""
# Issues: No methionine in sequence except for start
import numpy as np
import math
import scipy.stats as stats
import scipy.constants
import wc_kb
import wc_kb_gen
from Bio.Data import CodonTable
from Bio.Seq import Seq, Alphabet


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
    * translation_table (:obj: 'int'): The NCBI standard genetic code used
    * ncRNA_prop (:obj: 'float'): The proportion of non coding RNAs
    * rRNA_prop  (:obj: 'float'): The proportion of ribosomal RNAs
    * tRNA_prop (:obj: 'float'): The proportion of transfer RNAs
    * five_prime_len (:obj: 'int'): Average 5' UTR length for transcription units
    * three_prime_len (:obj: 'int'): Average 3' UTR length for transcription units
    * operon_prop (:obj: 'float'): Proportion of genes that should be in an operon (polycistronic mRNA)
    * operon_gen_num (:obj: 'int'): Average number of genes in an operon
    * mean_copy_number (:obj:`float`): mean copy number of each RNA
    * mean_half_life (:obj:`float`): mean half-life of each RNA in s
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

        # DOI: 10.1007/s10142-015-0433-4
        mean_coding_frac = options.get('mean_coding_frac', 0.88)
        assert(mean_coding_frac > 0 and mean_coding_frac < 1)
        options['mean_coding_frac'] = mean_coding_frac

        mean_num_genes = options.get('mean_num_genes', 4500)
        assert(mean_num_genes >= 1)
        options['mean_num_genes'] = mean_num_genes

        # default .53 doi: 10.1371/journal.pgen.1005613. p.11 2nd paragraph
        operon_prop = (options.get('operon_prop', 0.53))
        assert(operon_prop >= 0 and operon_prop <= 1)
        options['operon_prop'] = operon_prop

        mean_num_sRNA = options.get('mean_num_sRNA', 65)
        assert(mean_num_sRNA >= 0)
        options['mean_num_sRNA'] = mean_num_sRNA

        mean_num_rRNA = options.get('mean_num_rRNA', 25)
        assert(mean_num_rRNA >= 0)
        options['mean_num_rRNA'] = mean_num_rRNA

        mean_num_tRNA = options.get('mean_num_tRNA', 90)
        assert(mean_num_tRNA >= 0)
        options['mean_num_tRNA'] = mean_num_tRNA

        # Ensure the number of non mrna genes is less that the total number of monocistrocnic genes
        assert((mean_num_sRNA + mean_num_rRNA + mean_num_tRNA)
               < mean_num_genes * (1-operon_prop))

        # DOI: 10.1093/molbev/msk019
        # length in CODONS
        mean_gene_len = options.get('mean_gene_len', 308)
        assert(mean_gene_len >= 1)
        options['mean_gene_len'] = mean_gene_len

        translation_table = int(options.get('translation_table', 1))
        assert(translation_table in range(1, 32))
        options['translation_table'] = translation_table

        # default 150 DOI: 10.1016/S0968-0004(03)00002-1
        five_prime_len = int(options.get('five_prime_len', 13))
        assert(five_prime_len >= 12)
        options['five_prime_len'] = five_prime_len

        # default 500 DOI: 10.1016/S0968-0004(03)00002-1
        three_prime_len = int(options.get('three_prime_len', 5))
        assert(three_prime_len >= 0)
        options['three_prime_len'] = three_prime_len

        # default .53 doi: 10.1371/journal.pgen.1005613. p.11 2nd paragraph
        operon_prop = (options.get('operon_prop', 0.53))
        assert(operon_prop >= 0 and operon_prop <= 1)
        options['operon_prop'] = operon_prop

        min_operon_size = int(options.get('min_operon_size', 2))
        assert(min_operon_size >= 2)
        options['min_operon_size'] = min_operon_size

        max_operon_size = int(options.get('max_operon_size', 6))
        assert(max_operon_size >= min_operon_size)
        options['max_operon_size'] = max_operon_size

        # 1 or 4 bps doi: 10.1073/pnas.110147297
        operon_gen_spacing = (options.get('operon_gen_spacing', [1, 4]))
        shape = np.shape(operon_gen_spacing)
        assert(isinstance(shape, tuple))
        assert(len(shape) == 1)
        for i in shape:
            assert(isinstance(i, int))
            assert(i >= 1)
        options['operon_gen_spacing'] = operon_gen_spacing

        # DOI: 10.1038/ismej.2012.94
        mean_copy_number = options.get('mean_copy_number', 0.4)
        assert(mean_copy_number > 0)
        options['mean_copy_number'] = mean_copy_number

        # DOI: 10.1073/pnas.0308747101
        mean_half_life = options.get('mean_half_life', 2.1 * 60)
        assert(mean_half_life > 0)
        options['mean_half_life'] = mean_half_life

        # DOI: 10.1038/ismej.2012.94
        mean_copy_number = options.get('mean_copy_number', 0.4)
        assert(mean_copy_number > 0)
        options['mean_copy_number'] = mean_copy_number

        # DOI: 10.1073/pnas.0308747101
        mean_half_life = options.get('mean_half_life', 2.1 * 60)
        assert(mean_half_life > 0)
        options['mean_half_life'] = mean_half_life

        # DOI: 10.1038/ismej.2012.94
        mean_copy_number = options.get('mean_copy_number', 0.4)
        assert(mean_copy_number > 0)
        options['mean_copy_number'] = mean_copy_number

        # DOI: 10.1073/pnas.0308747101
        mean_half_life = options.get('mean_half_life', 2.1 * 60)
        assert(mean_half_life > 0)
        options['mean_half_life'] = mean_half_life

    def gen_components(self):
        np.random.seed()  # TODO match this with config
        self.gen_genome()
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
        mean_num_sRNA = options.get('mean_num_sRNA')
        mean_num_rRNA = options.get('mean_num_rRNA')
        mean_num_tRNA = options.get('mean_num_tRNA')
        min_operon_size = options.get('min_operon_size')
        max_operon_size = options.get('max_operon_size')
        operon_prop = options.get('operon_prop')
        three_prime_len = options.get('three_prime_len')
        five_prime_len = options.get('five_prime_len')
        operon_gen_spacing = options.get('operon_gen_spacing')
        cell = self.knowledge_base.cell

        # Define DNA Constants
        self.knowledge_base.translation_table = translation_table
        codon_table = CodonTable.unambiguous_dna_by_id[translation_table]
        START_CODONS = codon_table.start_codons
        STOP_CODONS = codon_table.stop_codons
        BASES = ['A', 'C', 'G', 'T']
        PROB_BASES = [(1 - mean_gc_frac) / 2, mean_gc_frac /
                      2, mean_gc_frac/2, (1-mean_gc_frac)/2]
        # The probability of each base being selected randomly

        # Create a chromosome n times
        for i_chr in range(num_chromosomes):

            # Determine numbers of each type of gene
            num_genes = self.rand(mean_num_genes/num_chromosomes)[0]

            num_polycistronic_genes = self.rand(
                operon_prop*num_genes, max=num_genes-1)[0]

            num_monocistronic_genes = num_genes-num_polycistronic_genes

            # Create a list of the size of each gene in BPS
            gene_lens = self.rand(mean_gene_len*3, count=num_genes)

            # Create a list of the size of each operon
            mean_operon_size = (min_operon_size+max_operon_size)/2
            operon_sizes = []

            while(np.sum(operon_sizes) < num_polycistronic_genes):
                operon_size_i = self.rand(
                    mean_operon_size, min=min_operon_size, max=max_operon_size)[0]
                if (np.sum(operon_sizes)+operon_size_i > num_polycistronic_genes):
                    operon_size_i = num_polycistronic_genes - \
                        np.sum(operon_sizes)
                operon_sizes.append(int(operon_size_i))

            # Determine num of operons and transcription units
            num_operons = len(operon_sizes)
            num_transcription_units = num_monocistronic_genes + num_operons

            # Begin gen sequence
            # total number of codons that code
            total_coding_len = np.sum(gene_lens)

            # avg number of codons that dont code
            mean_total_non_coding_len = self.rand(
                (1-mean_coding_frac)*total_coding_len)[0]

            # a list of non coding buffer between each transctiption unit

            inter_tu_lens = self.rand((mean_total_non_coding_len /
                                       num_transcription_units), count=num_transcription_units)

            total_non_coding_len = np.sum(inter_tu_lens)
            seq_len = 3*(total_coding_len + total_non_coding_len)

            seq_str = []
            # generates seq based on random codons (NOT start/stop codons)
            for i in range(0, seq_len, 3):
                codon_i = STOP_CODONS[0]

                while(codon_i in (STOP_CODONS or START_CODONS)):
                    codon_i = "".join(np.random.choice(
                        BASES, p=PROB_BASES, size=(3,)))

                seq_str.append(codon_i)

            seq_str = "".join(seq_str)

            # End Gen Seq

            chro = cell.species_types.get_or_create(
                id='chr_{}'.format(i_chr + 1), __type=wc_kb.DnaSpeciesType)
            chro.name = 'Chromosome {}'.format(i_chr + 1)
            chro.circular = chromosome_topology == 'circular'
            chro.double_stranded = True
            seq = Seq(seq_str, Alphabet.DNAAlphabet())
            chro.seq = seq

            tu_types = np.concatenate(
                (np.ones(num_monocistronic_genes), np.zeros(num_operons)),)
            np.random.shuffle(tu_types)

            loci = self.knowledge_base.cell.loci

            num_mRNA = 0

            while(num_mRNA < 1):
                num_sRNA = self.rand(mean_num_sRNA)[0]
                num_tRNA = self.rand(mean_num_tRNA)[0]
                num_rRNA = self.rand(mean_num_rRNA)[0]
                num_mRNA = num_monocistronic_genes-num_tRNA-num_rRNA-num_sRNA

            gene_types = [wc_kb.GeneType.mRna, wc_kb.GeneType.rRna,
                          wc_kb.GeneType.sRna, wc_kb.GeneType.tRna]

            a = np.full(num_mRNA, 0)
            b = np.full(num_rRNA, 1)
            c = np.full(num_sRNA, 2)
            d = np.full(num_tRNA, 3)
            gene_type_selector = np.concatenate((a, b, c, d))
            np.random.shuffle(gene_type_selector)
            gene_type_selector = iter(gene_type_selector)

            index = 13
            gene_len_iterator = enumerate(gene_lens)
            operon_size_iterator = iter(operon_sizes)

            def gen_gene(tu, allmRna=False):
                try:
                    i_gene, i_size = next(gene_len_iterator)
                except StopIteration:
                    return None, None

                gene = loci.get_or_create(id='gene_{}_{}'.format(
                    i_chr + 1, i_gene + 1), __type=wc_kb.GeneLocus)
                gene.name = 'Gene_{}_{}'.format(i_chr + 1, i_gene + 1)
                gene.start = index
                gene.end = index+i_size
                gene.polymer = chro

                gene.strand = tu.strand
                gene.transcription_units.append(tu)

                if (allmRna):  # True if called from operon construction
                    gene.type = wc_kb.GeneType.mRna
                else:  # If called from monocistrocnic construction
                    gene.type = gene_types[next(gene_type_selector)]

                if(gene.type == wc_kb.GeneType.mRna):  # If this gene is mRna
                    start_codon = np.random.choice(START_CODONS)
                    stop_codon = np.random.choice(STOP_CODONS)
                    seq_str = str(chro.seq)
                    seq_str = seq_str[: gene.start-1] + start_codon + \
                        seq_str[gene.start+2: gene.end-3] + \
                        stop_codon + seq_str[gene.end:]
                    chro.seq = Seq(seq_str, Alphabet.DNAAlphabet())

                index_i = index + i_size

                return gene, index_i

            for i_tu, tu_type in enumerate(tu_types):
                gene = None

                if tu_type:  # Is monocistronic

                    tu = loci.get_or_create(id='tu_{}_{}'.format(
                        i_chr + 1, i_tu + 1), __type=wc_kb.TranscriptionUnitLocus)
                    tu.polymer = chro
                    tu.name = 'Transcription unit {}-{}'.format(
                        i_chr + 1, i_tu + 1)
                    tu.strand = np.random.choice(
                        (wc_kb.PolymerStrand.positive, wc_kb.PolymerStrand.negative))

                    gene, index = gen_gene(tu)  # create a gene

                    # TODO fix this to ensure that this does not cause overlapping tus
                    tu.start = gene.start-self.rand(five_prime_len)[0]
                    three_prime = self.rand(three_prime_len)[0]
                    tu.end = gene.end+three_prime
                    index += three_prime
                else:
                    tu = loci.get_or_create(id='tu_{}_{}'.format(
                        i_chr + 1, i_tu + 1), __type=wc_kb.TranscriptionUnitLocus)
                    tu.polymer = chro
                    tu.name = 'Transcription unit {}-{}'.format(
                        i_chr + 1, i_tu + 1)
                    tu.strand = np.random.choice(
                        (wc_kb.PolymerStrand.positive, wc_kb.PolymerStrand.negative))

                    try:
                        for i in range(next(operon_size_iterator)):
                            gene, index = gen_gene(tu, True)  # create a gene
                            index += np.random.choice(operon_gen_spacing)
                    except StopIteration:
                        pass

                    # TODO fix this to ensure that this does not cause overlapping tus
                    tu.start = gene.start-self.rand(five_prime_len)[0]
                    three_prime = self.rand(three_prime_len)[0]
                    tu.end = gene.end+three_prime
                    index += three_prime

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

                tu = chromosome.loci[i]

                if type(tu) == wc_kb.TranscriptionUnitLocus:

                    # creates RnaSpeciesType for RNA sequence corresponding to gene

                    rna = self.knowledge_base.cell.species_types.get_or_create(
                        id='rna_{}'.format(tu.id), __type=wc_kb.RnaSpeciesType)
                    rna.name = 'rna {}'.format(tu.id)

                    if tu.genes[0].type == wc_kb.GeneType.mRna:
                        rna.type = wc_kb.RnaType.mRna
                    elif tu.genes[0].type == wc_kb.GeneType.rRna:
                        rna.type = wc_kb.RnaType.rRna
                    elif tu.genes[0].type == wc_kb.GeneType.tRna:
                        rna.type = wc_kb.RnaType.tRna
                    elif tu.genes[0].type == wc_kb.GeneType.sRna:
                        rna.type = wc_kb.RnaType.sRna

                    # print(rna.type)

                    rna.concentration = np.random.gamma(
                        1, mean_copy_number) / scipy.constants.Avogadro / mean_volume
                    rna.half_life = np.random.normal(
                        mean_half_life, np.sqrt(mean_half_life))

                    rna.transcription_units.append(tu)

                    # adds corresponding mRNA sequence to speciestypes list of kb.cell

                    if rna.type == wc_kb.RnaType.mRna:
                        for gene in tu.genes:
                            # creates ProteinSpeciesType object for corresponding protein sequence(s)
                            prot = self.knowledge_base.cell.species_types.get_or_create(
                                id='prot_{}'.format(gene.id), __type=wc_kb.ProteinSpeciesType)
                            prot.name = 'prot {}'.format(gene.id)
                            prot.gene = gene  # associates protein with GeneLocus object for corresponding gene
                            prot.rna = rna

    def rand(self, mean, count=1, min=0, max=np.inf):
        """ Generated 1 or more random normally distributed integer(s) between the minimum and maximum values,  with standard deviation equal  to the square root of the mean value.

        Args:
            mean (:obj:`float`): mean value
            count (:obj:`int`): number of random numbers to generate
            min   (:obj:'int'): the minimum value to generate
            max   (:obj:'int'): the maximum value to generate
        Returns:
           obj:`numpy.ndarray` of :obj:`int`: random normally distributed integer(s)
        """
        a = (min-mean)/np.sqrt(mean)
        b = (max - mean)/np.sqrt(mean)

        return np.int64(np.round(stats.truncnorm.rvs(a, b, loc=mean, scale=np.sqrt(mean), size=count)))
