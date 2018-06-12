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
import unittest
import wc_kb
import wc_lang


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['rand_wc_model_gen', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegexpMatches(context.Exception, 'usage: rand_wc_model_gen')

        with mock.patch('sys.argv', ['rand_wc_model_gen']):
            with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
                __main__.main()
                self.assertRegexpMatches(stdout.getvalue().strip(), 'usage: rand_wc_model_gen')
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

    def test_gen_model(self):
        kb_core_path = os.path.join(self.temp_dir, 'kb.xlsx')
        kb_seq_path = os.path.join(self.temp_dir, 'kb.fna')
        model_path = os.path.join(self.temp_dir, 'model.xlsx')
        with __main__.App(argv=['gen', kb_core_path, kb_seq_path, model_path]) as app:
            app.run()

        self.assertTrue(os.path.isfile(kb_core_path))
        self.assertTrue(os.path.isfile(kb_seq_path))
        self.assertTrue(os.path.isfile(model_path))
        kb = wc_kb.io.Reader().run(kb_core_path, kb_seq_path)
        model = wc_lang.io.Reader().run(model_path)

    def test_sim_model(self):
        kb_core_path = os.path.join(self.temp_dir, 'kb.xlsx')
        kb_seq_path = os.path.join(self.temp_dir, 'kb.fna')
        model_path = os.path.join(self.temp_dir, 'model.xlsx')
        with __main__.App(argv=['gen', kb_core_path, kb_seq_path, model_path]) as app:
            app.run()

        sim_results_path = os.path.join(self.temp_dir, 'sim_results')
        with __main__.App(argv=['sim', model_path, sim_results_path]) as app:
            app.run()

        # todo
