""" Test the generation of a random chromosome and accompanying mRNA and protein sequences
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""
import unittest
import numpy as np
from rand_wc_model_gen.kbgen.GenomeGenerator import GenomeGenerator

GEN_LEN = 20
INTER_LEN = 20
GEN_NUM = 20
TRANSLATION_TABLE = 1


class TestSynthetic(unittest.TestCase):
    def setUp(self):
        synthetic = Synthetic()

    def tearDown(self):
        pass

    def test_length(self):
        synthetic.generate(GEN_LEN, INTER_LEN, GEN_NUM, TRANSLATION_TABLE)

    def multiple_models(self):
        pass
