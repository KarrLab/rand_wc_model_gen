""" Generator for chromosomes, transcription units, and genes of random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from numpy import random
from Bio import Alphabet
from Bio import Seq
import numpy
import wc_kb
import wc_kb_gen


class ChromosomesGenesTusGenerator(wc_kb_gen.KbComponentGenerator):
    """ Generator for chromosomes, transcription units, and genes for random in silico organisms

    Options:

    * num_chromosomes (:obj:`int`): number of chromosomes
    * mean_gc_frac (:obj:`float`): fraction of chromosomes which are G or C
    * mean_num_genes (:obj:`float`): mean number of genes
    * mean_gene_len (:obj:`float`): mean length of a gene
    * mean_coding_frac (:obj:`float`): mean coding fraction of the genome
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        num_chromosomes = options.get('num_chromosomes', 1)
        assert(num_chromosomes >= 1 and int(
            num_chromosomes) == num_chromosomes)
        options['num_chromosomes'] = num_chromosomes

        chromosome_topology = options.get('chromosome_topology', 'circular')
        assert(chromosome_topology in ['circular', 'linear'])
        options['chromosome_topology'] = chromosome_topology

        mean_gc_frac = options.get('mean_gc_frac', 0.5)
        assert(mean_gc_frac >= 0 and mean_gc_frac <= 1)
        options['mean_gc_frac'] = mean_gc_frac

        mean_num_genes = options.get('mean_num_genes', 4377) # Escherichia coli K-12
        assert(mean_num_genes >= 1)
        options['mean_num_genes'] = mean_num_genes

        # DOI: 10.1093/molbev/msk019
        mean_gene_len = options.get('mean_gene_len', 924)
        assert(mean_gene_len >= 1)
        options['mean_gene_len'] = mean_gene_len

        # DOI: 10.1007/s10142-015-0433-4
        mean_coding_frac = options.get('mean_coding_frac', 0.88)
        assert(mean_coding_frac > 0 and mean_coding_frac < 1)
        options['mean_coding_frac'] = mean_coding_frac

    def gen_components(self):
        """ Construct knowledge base components """

        # get options
        options = self.options
        num_chromosomes = options.get('num_chromosomes')
        chromosome_topology = options.get('chromosome_topology')
        mean_gc_frac = options.get('mean_gc_frac')
        mean_num_genes = options.get('mean_num_genes')
        mean_gene_len = options.get('mean_gene_len')
        mean_coding_frac = options.get('mean_coding_frac')

        # generate chromosomes and genes
        cell = self.knowledge_base.cell
        for i_chr in range(num_chromosomes):
            num_genes = self.rand(mean_num_genes / num_chromosomes)[0]
            gene_lens = self.rand(mean_gene_len, count=num_genes)
            intergene_lens = self.rand(
                mean_gene_len / mean_coding_frac * (1 - mean_coding_frac), count=num_genes)

            seq_len = numpy.sum(gene_lens) + numpy.sum(intergene_lens)
            seq = Seq.Seq(''.join(random.choice(('A', 'C', 'G', 'T'), p=((1 - mean_gc_frac) / 2, mean_gc_frac /
                                                                         2, mean_gc_frac / 2, (1 - mean_gc_frac) / 2), size=(seq_len, ))), Alphabet.DNAAlphabet())

            chr = cell.species_types.get_or_create(
                id='chr_{}'.format(i_chr + 1), __type=wc_kb.core.DnaSpeciesType)
            chr.name = 'Chromosome {}'.format(i_chr + 1)
            chr.circular = chromosome_topology == 'circular'
            chr.double_stranded = True
            chr.seq = seq

            gene_starts = numpy.int64(numpy.cumsum(numpy.concatenate(([0], gene_lens[0:-1])) +
                                                   numpy.concatenate((numpy.round(intergene_lens[0:1] / 2), intergene_lens[1:]))))
            for i_gene in range(num_genes):
                tu = cell.loci.get_or_create(id='tu_{}_{}'.format(
                    i_chr + 1, i_gene + 1), __type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
                tu.polymer = chr
                tu.name = 'Transcription unit {}-{}'.format(
                    i_chr + 1, i_gene + 1)
                tu.start = gene_starts[i_gene]
                tu.end = gene_starts[i_gene] + gene_lens[i_gene] - 1
                tu.strand = random.choice(
                    (wc_kb.core.PolymerStrand.positive, wc_kb.core.PolymerStrand.negative))

                gene = cell.loci.get_or_create(id='gene_{}_{}'.format(
                    i_chr + 1, i_gene + 1), __type=wc_kb.prokaryote_schema.GeneLocus)
                gene.polymer = chr
                gene.transcription_units.append(tu)
                gene.name = 'Gene {}-{}'.format(i_chr + 1, i_gene + 1)
                gene.start = gene_starts[i_gene]
                gene.end = gene_starts[i_gene] + gene_lens[i_gene] - 1
                gene.type = wc_kb.core.GeneType.mRna
                gene.strand = tu.strand

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
