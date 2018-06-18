"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-06-06
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb
import wc_kb_gen
import math
import numpy as np
import random
from Bio.Seq import Seq


class GenomeGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates synthetic chromosome with randomized genes/intergenic regions. Creates RNA and protein objects corresponding to the genes on this chromosome. Associates the chromosome, RNAs, proteins
    with a knowledge base object (and its Cell attribute)

    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        #TODO: ASHWIN validate all new options

        options = self.options

        gen_len = int(options.get('gen_len', 300))  # for prokaryote (~924 bp)
        assert(gen_len > 0)
        options['gen_len'] = gen_len

        # for prokaryote (~100 bp)
        inter_len = int(options.get('inter_len', 30))
        assert(gen_len > 0)
        options['inter_len'] = inter_len

        # for E. coli (~4400); number of genes varies widely among prokaryotes
        gen_num = int(options.get('gen_num', 4400))
        assert(gen_num > 0)
        options['gen_num'] = gen_num

        translation_table = int(options.get('translation_table', 1))
        assert(translation_table in [1])
        options['translation_table'] = translation_table

    def gen_components(self):
        '''Construct knowledge base components'''   
         

        # create codon list
        # TODO BILAL enable use of translation table
        self.START_CODONS = ['ATG']  # start codon
        self.STOP_CODONS = ['TAG', 'TAA', 'TGA']  # stop codons

        self.knowledge_base.translation_table = translation_table

        # indexList of start/end positions of each gene, creates 'synthetic' chromosome
        self.indexList = self.gen_genome(
            gen_len, inter_len, gen_num)
        # creates RNA and protein objects corresponding to the genes on chromosome
        self.gen_rnas_proteins(gen_num, self.indexList)

    def gen_genome(self):
        """ Creates 'synthetic' chromsome with randomized genes/intergenic regions

        Args:
            gen_len (:obj:`int`): average gene length
            inter_len (:obj:`int`): average intergenic region length
            gen_num (:obj:`int`): number of genes on chromosome


        Returns:
            :obj:`list`: list of tuples of start and end positions of each gene on chromosome
        """

         # get options
        options = self.options
        gen_len = options.get('gen_len')
        gen_num = options.get('gen_num')
        translation_table = options.get('translation_table')
        mean_gc_frac = options.get('mean_gc_frac')
        mean_coding_frac = options.get('mean_coding_frac')
        num_chromosomes = options.get('num_chromosomes')



        #TODO BILAL Incorporate other generation function and account start/stop

        
        gene_dist = np.random.normal(gen_len, math.sqrt(gen_len), gen_num).tolist(
        )  # takes random samples out of Gaussian distribution with mean of average gene length
        gene_dist = [round(x) for x in gene_dist]
        # takes random samples out of Gaussian distribution with mean of average intergenic length
        inter_dist = np.random.normal(inter_len, math.sqrt(
            inter_len), gen_num).tolist()
        inter_dist = [round(x) for x in inter_dist]
        chromosome = wc_kb.DnaSpeciesType()
        seq = ''
        arr = ['A', 'G', 'C', 'T']
        indexList = []
        index = 1

        num_genes = self.rand(mean_num_genes / num_chromosomes)[0]
        gene_lens = self.rand(mean_gene_len, count=num_genes)
        intergene_lens = self.rand(mean_gene_len / mean_coding_frac * (1 - mean_coding_frac), count=num_genes)

        seq_len = numpy.sum(gene_lens) + numpy.sum(intergene_lens)
        seq = Seq.Seq(''.join(random.choice(('A', 'C', 'G', 'T'),
                                            p=((1 - mean_gc_frac) / 2, mean_gc_frac / 2, mean_gc_frac / 2, (1 - mean_gc_frac) / 2),
                                            size=(seq_len, ))),
                      Alphabet.DNAAlphabet())


        for i in range(2 * gen_num):
            if i % 2 == 0:  # if i is even, region is a gene
                gene = ''
                gen_length = gene_dist.pop()
                indexList.append((index, index + (gen_length * 3) - 1))
                gene += random.choice(self.START_CODONS)  # add start codon
                for k in range(gen_length - 2):
                    codon = self.STOP_CODONS[0]
                    # to make sure random codon is not stop codon (premature termination)
                    while codon in self.STOP_CODONS:
                        # create random codon
                        codon = random.choice(
                            arr) + random.choice(arr) + random.choice(arr)
                    gene += codon  # add randomly chosen codon
                gene += random.choice(self.STOP_CODONS)  # add stop codon
                index += gen_length * 3
                seq += gene

            else:  # if i is odd, region is intergenic
                inter_length = inter_dist.pop()
                index += inter_length * 3
                for k in range(inter_length):
                    # add randomly chosen base triple
                    seq += random.choice(arr) + \
                        random.choice(arr) + random.choice(arr)

        # associate the random chromosome sequence with the DnaSpeciesType object
        chromosome.seq = Seq(seq)
        # add chromosome to kb.cell speciestypes list

        self.knowledge_base.cell.species_types.append(chromosome)

        return indexList

    def gen_rnas_proteins(self, chromosomes):
        """ Creates RNA and protein objects corresponding to genes on chromosome

        Args:
            gen_num (:obj:`int`): number of genes on chromosome
            indexList (:obj: 'list'): list of tuples of start and end positions of each gene on chromosome

        """
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

        
