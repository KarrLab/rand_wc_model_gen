""" Tests of command line program

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import __main__
import abduct
import mock
import os
import rand_wc_model_gen
import shutil
import tempfile
import time
import unittest
import wc_kb
import wc_lang

@unittest.skip("broken_legacy")
class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.temp_dir, 'kb', 'core'))
        os.makedirs(os.path.join(self.temp_dir, 'kb', 'seq'))
        os.makedirs(os.path.join(self.temp_dir, 'model'))

        self.kb_core_path = os.path.join(self.temp_dir, 'kb', 'core', 'core.xlsx')
        self.kb_seq_path = os.path.join(self.temp_dir, 'kb', 'seq', 'seq.fna')
        self.model_path = os.path.join(self.temp_dir, 'model', 'model.xlsx')
        self.sim_results_path = os.path.join(self.temp_dir, 'sim_results')
        self.analysis_path = os.path.join(self.temp_dir, 'analysis')
        self.config_path = os.path.join(self.temp_dir, 'rand_wc_model_gen.cfg')

        # write configuration file
        with open(self.config_path, 'w') as file:
            file.write('[rand_wc_model_gen]\n')
            file.write('    [[kb_gen]]\n')
            file.write('        data_repo_metadata = false\n')
            file.write('        [[[component]]]\n')
            file.write('            [[[[GenomeGenerator]]]]\n')
            file.write('                num_chromosomes = 1\n')
            file.write('                mean_num_genes = 200.\n')
            file.write('                mean_gene_len = 10.\n')
            file.write('                seq_path = {}\n'.format(self.kb_seq_path))
            file.write('    [[kb]]\n')
            file.write('        [[[path]]]\n')
            file.write('            core = {}\n'.format(self.kb_core_path))
            file.write('            seq = {}\n'.format(self.kb_seq_path))
            file.write('    [[model_gen]]\n')
            file.write('        data_repo_metadata = false\n')
            file.write('    [[model]]\n')
            file.write('        path = {}\n'.format(self.model_path))
            file.write('    [[sim]]\n')
            file.write('        end_time = 0.1\n')
            file.write('        time_step = 0.01\n')
            file.write('    [[sim_results]]\n')
            file.write('        checkpoint_period = 0.01\n')
            file.write('        path = {}\n'.format(self.sim_results_path))
            file.write('    [[analysis]]\n')
            file.write('        path = {}\n'.format(self.analysis_path))

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['rand-wc-model-gen', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: rand-wc-model-gen')

        with mock.patch('sys.argv', ['rand-wc-model-gen']):
            with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
                __main__.main()
                self.assertRegex(stdout.getvalue().strip(), 'usage: rand-wc-model-gen')
                self.assertEqual(stderr.getvalue(), '')

    def test_get_version(self):
        with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
            with __main__.App(argv=['-v']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
            self.assertEqual(stdout.getvalue().strip(), rand_wc_model_gen.__version__)
            self.assertEqual(stderr.getvalue(), '')

        with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
            with __main__.App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
            self.assertEqual(stdout.getvalue().strip(), rand_wc_model_gen.__version__)
            self.assertEqual(stderr.getvalue(), '')

    @unittest.skipIf(os.getenv('CIRCLECI', '0') in ['1', 'true'], 'Too long for CircleCI')
    def test_generate(self):
        # verify that KB and model files haven't been created
        self.assertFalse(os.path.isfile(self.kb_core_path))
        self.assertFalse(os.path.isfile(self.kb_seq_path))
        self.assertFalse(os.path.isfile(self.model_path))

        # generate model
        with __main__.App(argv=['generate', '--config-path', self.config_path]) as app:
            app.run()

        self.assertTrue(os.path.isfile(self.kb_core_path))
        self.assertTrue(os.path.isfile(self.kb_seq_path))
        self.assertTrue(os.path.isfile(self.model_path))
        kb = wc_kb.io.Reader().run(self.kb_core_path, self.kb_seq_path)[wc_kb.KnowledgeBase][0]
        model = wc_lang.io.Reader().run(self.model_path)[wc_lang.Model][0]

    @unittest.skipIf(os.getenv('CIRCLECI', '0') in ['1', 'true'], 'Too long for CircleCI')
    def test_simulate(self):
        # generate model
        with __main__.App(argv=['generate', '--config-path', self.config_path]) as app:
            app.run()

        # simulate model
        with __main__.App(argv=['simulate', '--config-path', self.config_path]) as app:
            app.run()

        # assert
        self.assertTrue(os.path.isfile(os.path.join(app.results['sim_results_path'], 'run_results.h5')))

    @unittest.skipIf(os.getenv('CIRCLECI', '0') in ['1', 'true'], 'Too long for CircleCI')
    def test_analyze(self):
        # generate model
        with __main__.App(argv=['generate', '--config-path', self.config_path]) as app:
            app.run()

        # simulate model
        for i_sim in range(3):
            with __main__.App(argv=['simulate', '--config-path', self.config_path, '--seed', str(i_sim)]) as app:
                app.run()
            time.sleep(1.)  # todo: remove after results directory naming is fixed

        # analyze simulation results
        with __main__.App(argv=['analyze', '--config-path', self.config_path]) as app:
            app.run()
