""" Analysis of random whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-14
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_analysis
from .sim.rna import RnaSimulationAnalysis


class AnalysisRunner(wc_analysis.AnalysisRunner):
    DEFAULT_ANALYSES = (
        RnaSimulationAnalysis,
    )
