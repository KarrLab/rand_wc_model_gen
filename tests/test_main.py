""" Tests of command line program

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen import __main__
import abduct
import mock
import unittest
import rand_wc_model_gen


class CliTestCase(unittest.TestCase):

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
