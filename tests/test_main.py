""" Tests of command line program

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from random_wc_model_generator import __main__
import abduct
import mock
import unittest
import random_wc_model_generator


class CliTestCase(unittest.TestCase):

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['random_wc_model_generator', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegexpMatches(context.Exception, 'usage: random_wc_model_generator')

        with mock.patch('sys.argv', ['random_wc_model_generator']):
            with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
                __main__.main()
                self.assertRegexpMatches(stdout.getvalue().strip(), 'usage: random_wc_model_generator')
                self.assertEqual(stderr.getvalue(), '')

    def test_get_version(self):
        with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
            with __main__.App(argv=['-v']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
            self.assertEqual(stdout.getvalue().strip(), random_wc_model_generator.__version__)
            self.assertEqual(stderr.getvalue(), '')

        with abduct.captured(abduct.out(), abduct.err()) as (stdout, stderr):
            with __main__.App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
            self.assertEqual(stdout.getvalue().strip(), random_wc_model_generator.__version__)
            self.assertEqual(stderr.getvalue(), '')
