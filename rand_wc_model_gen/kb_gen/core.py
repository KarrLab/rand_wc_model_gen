""" Generator for KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .genome import GenomeGenerator
from .metabolites import MetabolitesGenerator
from .properties import PropertiesGenerator
from .observables import ObservablesGenerator
from numpy import random
import rand_wc_model_gen
import wc_kb
import wc_kb_gen


class KbGenerator(wc_kb_gen.KbGenerator):
    """ Generator for KBs for random in silico organisms

    * Circular chromosome

    Options:

    * id
    * version
    * component

        * GenomeGenerator
        * MetabolitesGenerator
        * PropertiesGenerator
    """

    DEFAULT_COMPONENT_GENERATORS = (
        PropertiesGenerator,
        GenomeGenerator,
        MetabolitesGenerator,
        ObservablesGenerator,
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
        options = self.options

        id = options.get('id', 'rand_wc_model')
        assert(isinstance(id, str) or id is None)
        options['id'] = id

        name = options.get('name', 'Random whole-cell model')
        assert(isinstance(name, str) or name is None)
        options['name'] = name

        version = options.get('version', rand_wc_model_gen.__version__)
        assert(isinstance(version, str) or version is None)
        options['version'] = version

        seed = options.get('seed', None)
        if seed is not None:
            try:
                int(seed)
            except ValueError:
                raise ValueError(
                    '`seed` option must be convertible to a 32 bit unsigned integer or `None`')
        options['seed'] = seed
