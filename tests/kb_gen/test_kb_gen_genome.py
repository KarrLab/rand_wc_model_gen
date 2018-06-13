""" Test the generation of a random chromosome and accompanying mRNA and protein sequences
:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-13
:Copyright: 2018, Karr Lab
:License: MIT
"""
import wc_kb
import unittest
import numpy as np
from rand_wc_model_gen.kbgen import genome

GEN_LEN = 20
INTER_LEN = 20
GEN_NUM = 20
TRANSLATION_TABLE = 1
gen = None


class TestSynthetic(unittest.TestCase):
    def setUp(self):
        kb = wc_kb.knowledge_base
       
        gen = genome.GenomeGenerator(kb)

    def tearDown(self):
        pass

    def test_run(self):
        gen.generate_genome(GEN_LEN, INTER_LEN, GEN_NUM, TRANSLATION_TABLE)
       

    def test_length(self):
        gen.generate_genome(GEN_LEN, INTER_LEN, GEN_NUM, TRANSLATION_TABLE)

    def multiple_models(self):
        pass
