"""
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Date: 2018-07-18
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb_gen


class ObservablesGenerator(wc_kb_gen.KbComponentGenerator):
    """
    Creates a set of  observables from the existing genes, proteins, and complexes. These observables are used in the random model generator as reaction rate modifiers.

    Options:
    * observables_list (:obj: 'list'): a list of the observables to be created
    """

    def clean_and_validate_options(self):
        """Apply default options and validate options"""
