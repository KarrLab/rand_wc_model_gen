"""
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
::Date: 2018-07-19
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb_gen


class CompartmentsGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates compartments for the knowledge base from the provided list of compartment ids and names.At least two compartments, with ids(names) 'c' (cytosol) and 'e' (extracellular), must be created.

    Options:
    * compartments (:obj:'list') a list of (:obj:'tuple') of type (:obj:'string',:obj:'string') that contain the id and name of the compartment.

    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options
        compartments = options.get(
            'compartments', [('c', 'cytosol'), ('e', 'extracellular space')])

        assert len(compartments) >= 2, "Must define at least 2 compartments"
        assert('c', 'cytosol') in compartments, "Could not find cytosol"
        assert (
            'e', 'extracellular space') in compartments, "Could not find extracellular space"

        options['compartments'] = compartments

    def gen_components(self):
        cell = self.knowledge_base.cell
        compartments = self.options.get('compartments')
        for compartment_id, compartment_name in compartments:
            compartment = cell.compartments.get_or_create(id=compartment_id)
            compartment.name = compartment_name
