"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-06-06
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb
import wc_kb_gen
import numpy
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
    * mean_gene_len (:obj:`float`): mean length of a gene
    * mean_coding_frac (:obj:`float`): mean coding fraction of the genome
    * translation_table (:obj: 'int'): The NCBI standard genetic code used
    * mean_num_ncRNA (:obj: 'int'): The number of non coding RNAs
    * mean_num_rRNA  (:obj: 'int'): The number of ribosomal RNAs
    * mean_num_tRNA (:obj: 'int'): The number of transfer RNAs
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

        mean_num_ncRNA = options.get('mean_num_ncRNA', 65)
        assert(mean_num_ncRNA >= 0)
        options['mean_num_ncRNA'] = mean_num_ncRNA

        mean_num_rRNA = options.get('mean_num_rRNA', 25)
        assert(mean_num_rRNA >= 0)
        options['mean_num_ncRNA'] = mean_num_rRNA

        mean_num_tRNA = options.get('mean_num_tRNA', 90)
        assert(mean_num_tRNA >= 0)
        options['mean_num_tRNA'] = mean_num_tRNA

        assert((mean_num_ncRNA + mean_num_rRNA + mean_num_tRNA) < mean_num_genes)

        # DOI: 10.1093/molbev/msk019
        mean_gene_len = options.get('mean_gene_len', 924)
        assert(mean_gene_len >= 1)
        options['mean_gene_len'] = mean_gene_len

        # DOI: 10.1007/s10142-015-0433-4
        mean_coding_frac = options.get('mean_coding_frac', 0.88)
        assert(mean_coding_frac > 0 and mean_coding_frac < 1)
        options['mean_coding_frac'] = mean_coding_frac

        translation_table = int(options.get('translation_table', 1))
        assert(translation_table in range(1, 32))
        options['translation_table'] = translation_table

    def gen_components(self):
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

        cell = self.knowledge_base.cell

        codon_table = self.knowledge_base.translation_table = CodonTable.unambiguous_dna_by_id[
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
            num_genes = self.rand(mean_num_genes / num_chromosomes)[0]
            gene_lens = self.rand(mean_num_genes, count=num_genes)

            intergene_lens = self.rand(
                mean_gene_len / mean_coding_frac * (1 - mean_coding_frac), count=num_genes)

            seq_len = numpy.sum(gene_lens) + numpy.sum(intergene_lens)

            seq_str = []

            for i in range(0, seq_len, 3):
                codon_i = STOP_CODONS[0]

                while(codon_i in STOP_CODONS):
                    codon_i = "".join(random.choice(
                        BASES, p=PROB_BASES, size=(3,)))

                seq_str.append(codon_i)

            seq_str = "".join(seq_str)
            seq = Seq(seq_str, Alphabet.DNAAlphabet)

            chr = cell.species_types.get_or_create(
                id='chr_{}'.format(i_chr + 1), __type=wc_kb.DnaSpeciesType)
            chr.name = 'Chromosome {}'.format(i_chr + 1)
            chr.circular = chromosome_topology == 'circular'
            chr.double_stranded = True
            chr.seq = seq

            gene_starts = numpy.int64(numpy.cumsum(numpy.concatenate(([0], gene_lens[0:-1])) +
                                                   numpy.concatenate((numpy.round(intergene_lens[0:1] / 2), intergene_lens[1:]))))

            # TODO: BILAL label gene loci and create objects. Determine how to pass this information to TU generation method
            self.knowledge_base.cell.species_types.append(chr)

    def gen_rnas_proteins(self, chromosomes):
        """ Creates RNA and protein objects corresponding to genes on chromosome

        Args:
            gen_num (:obj:`int`): number of genes on chromosome
            indexList (:obj: 'list'): list of tuples of start and end positions of each gene on chromosome

        """
<<<<<<< HEAD
        for chromosome in chromosomes:
            for locus in chromosome.loci:
                if type(locus) == wc_kb.TranscriptionUnitLocus():
                    # creates RnaSpeciesType for RNA sequence corresponding to gene
                    rna = wc_kb.RnaSpeciesType()
                    # GeneLocus object for gene sequence, attribute of ProteinSpeciesType object
                    tu = locus
                    
                    rna.type = tu.genes[0].type

                    rna.transcription_units.append(tu)

                    # adds corresponding mRNA sequence to speciestypes list of kb.cell
                    self.knowledge_base.cell.species_types.append(rna)
                    if rna.type == wc_kb.RnaType.mRna:
                        for gene in tu.genes:
                            # creates ProteinSpeciesType object for corresponding protein sequence(s)
                            prot = wc_kb.ProteinSpeciesType()

                            prot.cell = self.knowledge_base.cell
                            prot.cell.knowledge_base = self.knowledge_base
                            
                            prot.gene = gene  # associates protein with GeneLocus object for corresponding gene
                            prot.rna = rna

                            # adds ProteinSpeciesType object to kb.cell speciestypes list
                            self.knowledge_base.cell.species_types.append(prot)



    def make_tus(self, chromosomes):

        #TranscriptionUnitLocus: entire transcription unit
        #GeneLocus: part of transcription unit that is translated in the end

        #validate these options in the options method
        options = self.options
        five_prime_len = options.get('five_prime_len') #7 bp default (E. coli, wikipedia)
        three_prime_len = options.get('three_prime_len') #5 bp default guess
        operon_prop = options.get('operon_prop') #0.2 default guess 
        operon_gen_num = options.get('operon_gen_num') #3 genes default (https://academic.oup.com/gbe/article/5/11/2242/653613)
        
        for chromosome in chromosomes:
            seq = chromosome.seq
            i = 0
            while i < len(chromosome.loci): #PolymerSpeciesType also has loci list
                gene = chromosome.loci[i]
                if gene.type == wc_kb.GeneType.mRna:
                    #polycistronic mRNA (multiple GeneLocus objects per TranscriptionUnitLocus)

                    five_prime = round(np.random.normal(five_prime_len, math.sqrt(five_prime_len), 1).tolist()[0])
                    three_prime = round(np.random.normal(three_prime_len, math.sqrt(three_prime_len), 1).tolist()[0])

                    operon_prob = random.random()

                    if operon_prob <= operon_prop: #make an operon (polycistronic mRNA, put multiple genes in one TransUnitLocus)
                        operon_genes = round(np.random.normal(operon_gen_num, math.sqrt(operon_gen_num), 1).tolist()[0])
                        #add 3', 5' UTRs to the ends of the transcription unit (upstream of first gene, downstream of last gene)
                        tu = wc_kb.TranscriptionUnitLocus()
                        five_prime_start = gene.start - five_prime
                        if five_prime_start < 0:
                            five_prime_start = 0
                        tu.genes.append(gene)
                        tu.start = five_prime_start
                        tu.polymer = gene.polymer
                        for k in range(operon_genes-1):
                            i += 1
                            if i < len(chromosome.loci):
                                if (chromosome.loci[i]).type == wc_kb.GeneType.mRna:
                                    gene = chromosome.loci[i]
                                    tu.genes.append(gene)

                                else:
                                    break
                                
                            else:
                                break

                        three_prime_end = gene.end + three_prime
                        if three_prime_end >= len(seq):
                            three_prime_end = len(seq) - 1
                        tu.end = three_prime_end
                            

                    else: #make an individual transcription unit for the gene
                        five_prime_start = gene.start - five_prime
                        three_prime_end = gene.end + three_prime
                        if five_prime_start < 0:
                            five_prime_start = 0
                        if three_prime_end >= len(seq):
                            three_prime_end = len(seq) - 1
                        tu = wc_kb.TranscriptionUnitLocus()
                        tu.start = five_prime_start
                        tu.end = three_prime_end
                        tu.polymer = gene.polymer
                        tu.genes.append(gene)
                        chromosome.loci.append(tu)
                    
                else:
                    tu = wc_kb.TranscriptionUnitLocus()
                    tu.start = gene.start
                    tu.end = gene.end
                    tu.polymer = gene.polymer
                    tu.genes.append(gene)
                    chromosome.loci.append(tu)

                i += 1

        
=======
        # TODO ASHWIN work on TUs rather than genes

        chromosome = self.knowledge_base.cell.species_types[0]

        for i in range(gen_num):
            # creates RnaSpeciesType for RNA sequence corresponding to gene
            rna = wc_kb.RnaSpeciesType()
            # GeneLocus object for gene sequence, attribute of ProteinSpeciesType object
            gene = wc_kb.GeneLocus()
            gene.start = indexList[i][0]
            gene.end = indexList[i][1]
            gene.polymer = chromosome
            gene.cell = self.knowledge_base.cell
            # TranscriptionUnitLocus object - attribute of RnaSpeciesType object, associated with gene sequence
            transcription_locus = wc_kb.TranscriptionUnitLocus()
            transcription_locus.polymer = chromosome
            transcription_locus.start = gene.start
            transcription_locus.end = gene.end

            rna.transcription_units.append(transcription_locus)

            # adds corresponding mRNA sequence to speciestypes list of kb.cell
            self.knowledge_base.cell.species_types.append(rna)

            # creates ProteinSpeciesType object for corresponding protein sequence
            prot = wc_kb.ProteinSpeciesType()

            prot.cell = self.knowledge_base.cell
            prot.cell.knowledge_base = self.knowledge_base

            prot.gene = gene  # associates protein with GeneLocus object for corresponding gene

            # adds ProteinSpeciesType object to kb.cell speciestypes list
            self.knowledge_base.cell.species_types.append(prot)

    def rand(self, mean, count=1):
        """ Generated 1 or more random normally distributed integer(s) with standard deviation equal
        to the square root of the mean value.

        Args:
            mean (:obj:`float`): mean value
            count (:obj:`int`): number of random numbers to generate

        Returns:
            :obj:`int` or :obj:`numpy.ndarray` of :obj:`int`: random normally distributed integer(s)
        """
        return numpy.int64(numpy.round(random.normal(mean, numpy.sqrt(mean), (count, ))))
>>>>>>> 7bd9bd54edc9f6331818f2e2d44906abb16ab193
