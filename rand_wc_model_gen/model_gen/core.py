""" Generator for models based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .transcription import TranscriptionSubmodelGenerator
from .rna_degradation import RnaDegradationSubmodelGenerator
import wc_model_gen


class ModelGenerator(wc_model_gen.ModelGenerator):
    """ Generator for models based on KBs for random in silico organisms
    """

    DEFAULT_COMPONENT_GENERATORS = (
        TranscriptionSubmodelGenerator,
        #RnaDegradationSubmodelGenerator,
    )
