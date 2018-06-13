""" Generator for metabolites of random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import numpy
import wc_kb
import wc_kb_gen


class MetabolitesGenerator(wc_kb_gen.KbComponentGenerator):
    """ Generator for metabolites for random in silico organisms
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        mean_ntp_conc = float(options.get('mean_ntp_conc', 1.5e-3))  # DOI: 10.1038/srep06522
        assert(mean_ntp_conc > 0)
        options['mean_ntp_conc'] = mean_ntp_conc

        mean_h2o_conc = float(options.get('mean_h2o_conc', 55.))  # DOI: 10.1073/pnas.0308747101
        assert(mean_h2o_conc > 0)
        options['mean_h2o_conc'] = mean_h2o_conc

        mean_ph = float(options.get('mean_ph', 7.5))  # DOI: 10.1128/JB.00615-07
        assert(mean_ph > 0)
        options['mean_ph'] = mean_ph

    def gen_components(self):
        """ Construct knowledge base components """
        cell = self.knowledge_base.cell

        # get options
        options = self.options
        mean_ntp_conc = options.get('mean_ntp_conc')
        mean_h2o_conc = options.get('mean_h2o_conc')
        mean_ph = options.get('mean_ph')

        # generate metabolites
        cell.species_types.get_or_create(__type=wc_kb.MetaboliteSpeciesType, id='atp', name='ATP', concentration=mean_ntp_conc)
        cell.species_types.get_or_create(__type=wc_kb.MetaboliteSpeciesType, id='ctp', name='CTP', concentration=mean_ntp_conc)
        cell.species_types.get_or_create(__type=wc_kb.MetaboliteSpeciesType, id='gtp', name='GTP', concentration=mean_ntp_conc)
        cell.species_types.get_or_create(__type=wc_kb.MetaboliteSpeciesType, id='utp', name='UTP', concentration=mean_ntp_conc)
        cell.species_types.get_or_create(__type=wc_kb.MetaboliteSpeciesType, id='ppi', name='PPI', concentration=mean_ntp_conc)
        cell.species_types.get_or_create(__type=wc_kb.MetaboliteSpeciesType, id='h2o', name='H2O', concentration=mean_h2o_conc)
        cell.species_types.get_or_create(__type=wc_kb.MetaboliteSpeciesType, id='h', name='H', concentration=numpy.power(10, -mean_ph))
