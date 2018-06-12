""" Generator for chromosomes and genes of random in silico organisms

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


class ChromosomesGenesGenerator(wc_kb_gen.KbComponentGenerator):
    """ Generator for chromosomes and genes for random in silico organisms

    Options:

    * num_chromosomes (:obj:`int`): number of chromosomes
    * avg_gc_frac (:obj:`float`): fraction of chromosomes which are G or C
    * avg_num_genes (:obj:`float`): average number of genes
    * avg_gene_len (:obj:`float`): average length of a gene
    * avg_coding_frac (:obj:`float`): average coding fraction of the genome
    """

    def get_data(self):
        """ Get data for knowledge base components """
        pass  # pragma: no cover

    def process_data(self):
        """ Process data for knowledge base components """
        pass  # pragma: no cover

    def gen_components(self):
        """ Construct knowledge base components """

        # get and validate options
        num_chromosomes = self.options.get('num_chromosomes', 1)
        assert(num_chromosomes >= 1 and int(num_chromosomes) == num_chromosomes)

        avg_gc_frac = self.options.get('avg_gc_frac', 0.5)
        assert(avg_gc_frac >= 0 and avg_gc_frac <= 1)

        avg_num_genes = self.options.get('avg_num_genes', 2000)
        assert(avg_num_genes >= 1)

        avg_gene_len = self.options.get('avg_gene_len', 924)  # DOI: 10.1093/molbev/msk019
        assert(avg_gene_len >= 1)

        avg_coding_frac = self.options.get('avg_coding_frac', 0.88)  # DOI: 10.1007/s10142-015-0433-4
        assert(avg_coding_frac > 0 and avg_coding_frac < 1)

        # generate chromosomes and genes
        cell = self.knowledge_base.cell
        for i_chr in range(num_chromosomes):
            num_genes = self.rand(avg_num_genes / num_chromosomes)[0]
            gene_lens = self.rand(avg_gene_len, count=num_genes)
            intergene_lens = self.rand(avg_gene_len / avg_coding_frac * (1 - avg_coding_frac), count=num_genes)

            seq_len = numpy.sum(gene_lens) + numpy.sum(intergene_lens)
            seq = Seq.Seq(''.join(random.choice(('A', 'C', 'G', 'T'),
                                                p=((1 - avg_gc_frac) / 2, avg_gc_frac / 2, avg_gc_frac / 2, (1 - avg_gc_frac) / 2),
                                                size=(seq_len, ))),
                          Alphabet.DNAAlphabet())

            chr = wc_kb.DnaSpeciesType(
                cell=cell,
                id='chr_{}'.format(i_chr + 1),
                name='Chromosome {}'.format(i_chr + 1),
                circular=True,
                double_stranded=True,
                seq=seq)
            
            gene_starts = numpy.int64(numpy.cumsum(numpy.concatenate(([0], gene_lens[0:-1])) +
                                       numpy.concatenate((numpy.round(intergene_lens[0:1] / 2), intergene_lens[1:]))))
            for i_gene in range(num_genes):
                wc_kb.GeneLocus(
                    cell=cell,
                    polymer=chr,
                    id='gene_{}_{}'.format(i_chr + 1, i_gene + 1),
                    name='Gene {}-{}'.format(i_chr + 1, i_gene + 1),
                    start=gene_starts[i_gene],
                    end=gene_starts[i_gene] + gene_lens[i_gene] - 1
                )

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
