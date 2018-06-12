""" Generator for KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .chrs_genes import ChromosomesGenesGenerator
from .properties import PropertiesGenerator
import wc_kb_gen


class KbGenerator(wc_kb_gen.KbGenerator):
    """ Generator for KBs for random in silico organisms

    * Circular chromosome
    """

    DEFAULT_COMPONENT_GENERATORS = (
        PropertiesGenerator,
        ChromosomesGenesGenerator,        
    )
