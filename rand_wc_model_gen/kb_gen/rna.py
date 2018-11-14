""" Generator for RNA of random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from numpy import random
import numpy
import scipy.constants
import wc_kb
import wc_kb_gen


class RnaGenerator(wc_kb_gen.KbComponentGenerator):
    """ Generator for RNA for random in silico organisms

    Options:

    * mean_copy_number (:obj:`float`): mean copy number of each RNA
    * mean_half_life (:obj:`float`): mean half-life of each RNA in s
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        mean_copy_number = options.get('mean_copy_number', 0.4)  # DOI: 10.1038/ismej.2012.94
        assert(mean_copy_number > 0)
        options['mean_copy_number'] = mean_copy_number

        mean_half_life = options.get('mean_half_life', 2.1 * 60)  # DOI: 10.1073/pnas.0308747101
        assert(mean_half_life > 0)
        options['mean_half_life'] = mean_half_life

    def gen_components(self):
        """ Construct knowledge base components """
        cell = self.knowledge_base.cell

        # get options
        options = self.options
        mean_copy_number = options.get('mean_copy_number')
        mean_half_life = options.get('mean_half_life')
        mean_volume = cell.properties.get_one(id='mean_volume').value

        # generate RNA
        cytosol = cell.compartments.get_one(id='c')
        tus = cell.loci.get(__type=wc_kb.prokaryote_schema.TranscriptionUnitLocus)
        for tu in tus:
            rna = cell.species_types.get_or_create(id=tu.id.replace('tu_', 'rna_'), __type=wc_kb.prokaryote_schema.RnaSpeciesType)
            rna.transcription_units = [tu]
            rna.name = tu.name.replace('Transcription unit', 'RNA')
            rna.type = wc_kb.core.RnaType[tu.genes[0].type.name]
            rna.species.get_or_create(compartment=cytosol).concentration = wc_kb.core.Concentration(
                cell=cell,
                value=random.gamma(1, mean_copy_number) / scipy.constants.Avogadro / mean_volume)
            rna.half_life = random.normal(mean_half_life, numpy.sqrt(mean_half_life))
