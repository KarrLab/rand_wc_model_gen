""" Generator for metabolism submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen


class MetabolismSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for metabolism submodel """

    def gen_compartments(self):
        cell = self.knowledge_base.cell
        model = self.model

        cyt = model.compartments.get_or_create(id='c')
        cyt.name = 'cytosol'
        cyt.initial_volume = cell.properties.get_one(id='mean_volume').value

        ext = model.compartments.get_or_create(id='e')
        ext.name = 'extracellular space'
        ext.initial_volume = 1. / cell.properties.get_one(id='mean_cell_density').value

    def gen_parameters(self):
        cell = self.knowledge_base.cell
        model = self.model
        param = model.parameters.get_or_create(id='fractionDryWeight')
        param.submodels.append(self.submodel)
        param.value = cell.properties.get_one(id='mean_fraction_dry_weight').value
        param.units = 'dimensionless'

    def gen_species(self):
        pass

    def gen_reactions(self):
        pass

    def gen_rate_laws(self):
        pass
