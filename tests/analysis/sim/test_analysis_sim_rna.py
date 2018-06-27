""" Tests of analysis of knowledge bases, models, and simulation results

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-14
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import analysis
from rand_wc_model_gen import kb_gen
from rand_wc_model_gen import model_gen
from wc_utils.util import rand
import os
import shutil
import tempfile
import time
import unittest
import wc_sim.multialgorithm.simulation


class RnaSimulationAnalysisTestCase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test(self):
        # generate kb
        kb = kb_gen.KbGenerator(options={
            'component': {
                'GenomeGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 100.,
                    'mean_gene_len': 100.,
                    'mean_num_rRNA': 1,
                    'mean_num_tRNA': 1,
                    'mean_num_sRNA': 1,
                },
            },
        }).run()

        # generate model
        model = model_gen.ModelGenerator(kb).run()

        # run simulations
        sim_results_path = os.path.join(self.temp_dir, 'sim_results')
        for i_sim in range(3):
            sim = wc_sim.multialgorithm.simulation.Simulation(model)
            sim.run(end_time=10.,
                    seed=i_sim,
                    checkpoint_period=1.,
                    results_dir=sim_results_path)
            # todo: remove after results directory naming is fixed
            time.sleep(1.0)

        # analyze simulation results
        analysis_results_path = os.path.join(self.temp_dir, 'analysis')
        analysis.sim.rna.RnaSimulationAnalysis(
            kb, model, sim_results_path, out_path=analysis_results_path).run()

        self.assertTrue(os.path.isfile(os.path.join(
            analysis_results_path, 'Individual RNA (simulation 1).pdf')))
        self.assertTrue(os.path.isfile(os.path.join(
            analysis_results_path, 'Total RNA (all simulations).pdf')))
