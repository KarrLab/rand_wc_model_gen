""" Generator for KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .chrs_genes import ChromosomesGenesGenerator
from .properties import PropertiesGenerator
from .rna import RnaGenerator
from numpy import random
import wc_kb
import wc_kb_gen


class KbGenerator(wc_kb_gen.KbGenerator):
    """ Generator for KBs for random in silico organisms

    * Circular chromosome

    Options:

    * id
    * version
    * component

        * PropertiesGenerator
        * ChromosomesGenesGenerator
        * RnaGenerator
    """

    DEFAULT_COMPONENT_GENERATORS = (
        PropertiesGenerator,
        ChromosomesGenesGenerator,
        RnaGenerator,
    )

    def run(self):
        """ Generate a knowledge base of experimental data for a whole-cell model

        Returns:
            :obj:`wc_kb.KnowledgeBase`: knowledge base
        """
        self.clean_and_validate_options()
        random.seed(self.options.get('seed'))
        return super(KbGenerator, self).run()

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        super(KbGenerator, self).clean_and_validate_options()

        options = self.options

        seed = options.get('seed', None)
        if seed is not None:
            try:
                int(seed)
            except ValueError:
                raise ValueError('`seed` option must be convertible to a 32 bit unsigned integer or `None`')
        options['seed'] = seed
