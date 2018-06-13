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

    * mean_volume (:obj:`float`): mean volume in L
    * mean_doubling_time (:obj:`float`): mean doubling time in s
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        mean_volume = options.get('mean_volume', 1e-15)
        assert(mean_volume > 0)
        options['mean_volume'] = mean_volume

        mean_fraction_dry_weight = options.get('mean_fraction_dry_weight', 0.3)
        assert(mean_fraction_dry_weight > 0)
        options['mean_fraction_dry_weight'] = mean_fraction_dry_weight

        mean_doubling_time = options.get('mean_doubling_time', 30 * 60)
        assert(mean_doubling_time > 0)
        options['mean_doubling_time'] = mean_doubling_time

        mean_cell_density = options.get('mean_cell_density', 1e6)
        assert(mean_cell_density > 0)
        options['mean_cell_density'] = mean_cell_density

    def gen_components(self):
        """ Construct knowledge base components """

        # get options
        options = self.options

        # generate properties
        cell = self.knowledge_base.cell

        prop = cell.properties.get_or_create(id='mean_volume')
        prop.value = options.get('mean_volume')
        prop.units = 'L'

        prop = cell.properties.get_or_create(id='mean_fraction_dry_weight')
        prop.value = options.get('mean_fraction_dry_weight')
        prop.units = 'dimensionless'

        prop = cell.properties.get_or_create(id='mean_doubling_time')
        prop.value = options.get('mean_doubling_time')
        prop.units = 's'

        prop = cell.properties.get_or_create(id='mean_cell_density')
        prop.value = options.get('mean_cell_density')
        prop.units = 'cells L^{-1}'
