""" Tests of analysis of knowledge bases, models, and simulation results

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-14
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import analysis
from rand_wc_model_gen import kb_gen
import os
import shutil
import tempfile
import time
import unittest
import wc_model_gen.prokaryote
import wc_sim.multialgorithm.simulation

@unittest.skip("broken_legacy")
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
                    'mean_num_genes': 200.,
                    'mean_gene_len': 10.,
                    'seq_path': os.path.join(self.temp_dir, 'kb_seq.fna')
                },
            },
        }).run()

        # generate model
        model = wc_model_gen.prokaryote.ProkaryoteModelGenerator(kb).run()

        # run simulations
        sim_results_path = os.path.join(self.temp_dir, 'sim_results')
        for i_sim in range(3):
            sim = wc_sim.multialgorithm.simulation.Simulation(model)
            sim.run(end_time=1e-1,
                    time_step=1e-2,
                    seed=i_sim,
                    checkpoint_period=1e-2,
                    results_dir=sim_results_path)
            # todo: remove after results directory naming is fixed
            time.sleep(1.0)

        # analyze simulation results
        analysis_results_path = os.path.join(self.temp_dir, 'analysis')
        analysis.sim.rna.RnaSimulationAnalysis(
            sim_results_path, knowledge_base=kb, model=model, out_path=analysis_results_path).run()

        self.assertTrue(os.path.isfile(os.path.join(
            analysis_results_path, 'Individual RNA (simulation 1).pdf')))
        self.assertTrue(os.path.isfile(os.path.join(
            analysis_results_path, 'Total RNA (all simulations).pdf')))
