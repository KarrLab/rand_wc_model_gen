""" Generator for properties

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import wc_kb_gen


class PropertiesGenerator(wc_kb_gen.KbComponentGenerator):
    """ Generator for other properties for random in silico organisms

    Options:

    * mean_doubling_time (:obj:`float`): mean doubling time in s
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        mean_doubling_time = options.get('mean_doubling_time', 30 * 60)
        assert(mean_doubling_time > 0)
        options['mean_doubling_time'] = mean_doubling_time

    def gen_components(self):
        """ Construct knowledge base components """

        # get options
        options = self.options
        mean_doubling_time = options.get('mean_doubling_time')

        # generate properties
        cell = self.knowledge_base.cell
        prop = cell.properties.get_or_create(id='mean_doubling_time')
        prop.value = mean_doubling_time
        prop.units = 's'