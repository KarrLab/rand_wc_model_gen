""" Tests of analysis of knowledge bases, models, and simulation results

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-14
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import analysis
import unittest

@unittest.skip("broken_legacy")
class AnalysisRunnerTestCase(unittest.TestCase):
    def test(self):
        analysis.AnalysisRunner(None, None, None, analyses=())
