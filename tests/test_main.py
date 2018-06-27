""" Tests of command line program

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import __main__
from wc_utils.util import rand
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


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.kb_core_path = os.path.join(
            self.temp_dir, 'kb', 'core', 'core.xlsx')
        self.kb_seq_path = os.path.join(self.temp_dir, 'kb', 'seq', 'seq.fna')
        self.model_path = os.path.join(self.temp_dir, 'model', 'model.xlsx')
        self.sim_results_path = os.path.join(self.temp_dir, 'sim_results')
        self.analysis_path = os.path.join(self.temp_dir, 'analysis')
        self.config_path = os.path.join(self.temp_dir, 'rand_wc_model_gen.cfg')

        # write configuration file
        with open(self.config_path, 'w') as file:
            file.write('[rand_wc_model_gen]\n')
            file.write('    [[kb_gen]]\n')
            file.write('        set_repo_metadata_from_path = false\n')
            file.write('        [[[component]]]\n')
            file.write('            [[[[GenomeGenerator]]]]\n')
            file.write('                num_chromosomes = 1\n')
            file.write('                mean_num_rRNA = 5\n')
            file.write('                mean_num_tRNA = 5\n')
            file.write('                mean_num_sRNA = 5\n')
            file.write('                mean_num_genes = 100.\n')
            file.write('                mean_gene_len = 10.\n')
            file.write('    [[kb]]\n')
            file.write('        [[[path]]]\n')
            file.write('            core = {}\n'.format(self.kb_core_path))
            file.write('            seq = {}\n'.format(self.kb_seq_path))
            file.write('    [[model_gen]]\n')
            file.write('        set_repo_metadata_from_path = false\n')
            file.write('    [[model]]\n')
            file.write('        path = {}\n'.format(self.model_path))
            file.write('    [[sim]]\n')
            file.write('        end_time = 5.\n')
            file.write('    [[sim_results]]\n')
            file.write('        checkpoint_period = 1.\n')
            file.write('        path = {}\n'.format(self.sim_results_path))
            file.write('    [[analysis]]\n')
            file.write('        path = {}\n'.format(self.analysis_path))

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['rand_wc_model_gen', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegexpMatches(
                    context.Exception, 'usage: rand_wc_model_gen')

        with mock.patch('sys.argv', ['rand_wc_model_gen']):
            with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
                __main__.main()
                self.assertRegexpMatches(
                    stdout.getvalue().strip(), 'usage: rand_wc_model_gen')
                self.assertEqual(stderr.getvalue(), '')

    def test_get_version(self):
        with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
            with __main__.App(argv=['-v']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
            self.assertEqual(stdout.getvalue().strip(),
                             rand_wc_model_gen.__version__)
            self.assertEqual(stderr.getvalue(), '')

        with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
            with __main__.App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
            self.assertEqual(stdout.getvalue().strip(),
                             rand_wc_model_gen.__version__)
            self.assertEqual(stderr.getvalue(), '')

    def test_generate(self):
        # generate model
        with __main__.App(argv=['generate', '--config-path', self.config_path]) as app:
            app.run()

        self.assertTrue(os.path.isfile(self.kb_core_path))
        self.assertTrue(os.path.isfile(self.kb_seq_path))
        self.assertTrue(os.path.isfile(self.model_path))
        kb = wc_kb.io.Reader().run(self.kb_core_path, self.kb_seq_path)
        model = wc_lang.io.Reader().run(self.model_path)

    def test_simulate(self):
        # generate model
        with __main__.App(argv=['generate', '--config-path', self.config_path]) as app:
            app.run()

        # simulate model
        with __main__.App(argv=['simulate', '--config-path', self.config_path]) as app:
            app.run()

        # assert
        self.assertTrue(os.path.isfile(os.path.join(
            app.results['sim_results_path'], 'run_results.h5')))

    def test_analyze(self):
        # generate model
        with __main__.App(argv=['generate', '--config-path', self.config_path]) as app:
            app.run()

        # simulate model
        for i_sim in range(3):
            with __main__.App(argv=['simulate', '--config-path', self.config_path, '--seed', str(i_sim)]) as app:
                app.run()
            # todo: remove after results directory naming is fixed
            time.sleep(1.)

        # analyze simulation results
        with __main__.App(argv=['analyze', '--config-path', self.config_path]) as app:
            app.run()
