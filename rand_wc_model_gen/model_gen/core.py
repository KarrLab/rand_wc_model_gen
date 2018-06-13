""" Generator for models based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .metabolism import MetabolismSubmodelGenerator
from .rna_degradation import RnaDegradationSubmodelGenerator
from .transcription import TranscriptionSubmodelGenerator
import wc_model_gen


class ModelGenerator(wc_model_gen.ModelGenerator):
    """ Generator for models based on KBs for random in silico organisms
    """

    DEFAULT_COMPONENT_GENERATORS = (
        MetabolismSubmodelGenerator,
        TranscriptionSubmodelGenerator,
        #RnaDegradationSubmodelGenerator,
    )
