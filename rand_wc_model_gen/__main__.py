""" Command line programs for generating random whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
import rand_wc_model_gen
import rand_wc_model_gen.config
import rand_wc_model_gen.kb_gen
import rand_wc_model_gen.model_gen
import wc_kb.io
import wc_lang.io
import wc_sim.multialgorithm


class BaseController(CementBaseController):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "Random whole-cell model generator"
        arguments = [
            (['-v', '--version'], dict(action='version', version=rand_wc_model_gen.__version__)),
        ]

    @expose(hide=True)
    def default(self):
        self.app.args.print_help()


class GenController(CementBaseController):
    """ Generate random whole-cell model """

    class Meta:
        label = 'gen'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Generate random whole-cell model"
        arguments = [
            (['kb_core_path'], dict(type=str, help='Path to save knowledge base core (.csv, .tsv, .xlsx)')),
            (['kb_seq_path'], dict(type=str, help='Path to save knowledge base sequence (.fna)')),
            (['model_path'], dict(type=str, help='Path to save model (.csv, .tsv, .xlsx)')),
            (['--config-path'], dict(type=str, default=None, help='Path to configuration file')),
            (['--set-repo-metadata-from-path'], dict(type=str, default=None, help='Path to configuration file')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        config = rand_wc_model_gen.config.get_config(extra_path=args.config_path)['rand_wc_model_gen']
        kb = rand_wc_model_gen.kb_gen.KbGenerator(options=config['kb_gen']).run()
        model = rand_wc_model_gen.model_gen.ModelGenerator(kb, options=config['model_gen']).run()
        wc_kb.io.Writer().run(kb, args.kb_core_path, args.kb_seq_path, set_repo_metadata_from_path=args.set_repo_metadata_from_path)
        wc_lang.io.Writer().run(model, args.model_path, set_repo_metadata_from_path=args.set_repo_metadata_from_path)


class SimController(CementBaseController):
    """ Simulate random whole-cell model """

    class Meta:
        label = 'sim'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Simulate random whole-cell model"
        arguments = [
            (['model_path'], dict(type=str, help='Path to model (.csv, .tsv, .xlsx)')),
            (['sim_results_path'], dict(type=str, help='Path to save simulation results in HDF5 format')),
            (['--len'], dict(type=float, help='Length of time to simulation in s')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        model = wc_lang.io.Reader().run(args.model_path)
        # wc_sim.multialgorithm(model, args.results_path) # todo


class App(CementApp):
    """ Command line application """
    class Meta:
        label = 'rand_wc_model_gen'
        base_controller = 'base'
        handlers = [
            BaseController,
            GenController,
            SimController,
        ]


def main():
    with App() as app:
        app.run()
