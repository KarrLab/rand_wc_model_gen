""" Command line programs for generating random whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
import glob
import os
import rand_wc_model_gen
import rand_wc_model_gen.analysis
import rand_wc_model_gen.config
import rand_wc_model_gen.kb_gen
import rand_wc_model_gen.model_gen
import wc_kb.io
import wc_lang.io
import wc_sim.multialgorithm.run_results
import wc_sim.multialgorithm.simulation


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


class GenerateController(CementBaseController):
    """ Generate a random whole-cell knowledge base and a random whole-cell model """

    class Meta:
        label = 'generate'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Generate a random whole-cell knowledge base and a random whole-cell model"
        arguments = [
            (['--config-path'], dict(type=str,
                                     default=None,
                                     help='Path to configuration file')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        config = rand_wc_model_gen.config.get_config(extra_path=args.config_path)['rand_wc_model_gen']
        kb = rand_wc_model_gen.kb_gen.KbGenerator(options=config['kb_gen']).run()
        model = rand_wc_model_gen.model_gen.ModelGenerator(kb, options=config['model_gen']).run()

        if not os.path.isdir(os.path.dirname(config['kb']['path']['core'])):
            os.makedirs(os.path.dirname(config['kb']['path']['core']))
        if not os.path.isdir(os.path.dirname(config['kb']['path']['seq'])):
            os.makedirs(os.path.dirname(config['kb']['path']['seq']))
        wc_kb.io.Writer().run(kb,
                              config['kb']['path']['core'], config['kb']['path']['seq'],
                              set_repo_metadata_from_path=config['kb_gen']['set_repo_metadata_from_path'])

        if not os.path.isdir(os.path.dirname(config['model']['path'])):
            os.makedirs(os.path.dirname(config['model']['path']))
        wc_lang.io.Writer().run(model,
                                config['model']['path'],
                                set_repo_metadata_from_path=config['model_gen']['set_repo_metadata_from_path'])


class SimulateController(CementBaseController):
    """ Simulate a random whole-cell model """

    class Meta:
        label = 'simulate'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Simulate a random whole-cell model"
        arguments = [
            (['--config-path'], dict(type=str, default=None, help='Path to configuration file')),
            (['--seed'], dict(type=int, default=None, help='Random number generator seed')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        config_path = args.config_path
        seed = args.seed
        if seed is None:
            seed = config['sim']['seed']

        config = rand_wc_model_gen.config.get_config(extra_path=config_path)['rand_wc_model_gen']
        model = wc_lang.io.Reader().run(config['model']['path'])
        simulation = wc_sim.multialgorithm.simulation.Simulation(model)
        num_events, sim_results_path = simulation.run(end_time=config['sim']['end_time'],
                                                      seed=seed,
                                                      checkpoint_period=config['sim_results']['checkpoint_period'],
                                                      results_dir=config['sim_results']['path'])
        self.app.results = {
            'num_events': num_events,
            'sim_results_path': sim_results_path,
        }


class AnalyzeController(CementBaseController):
    """ Analyze a random whole-cell knowledge base, model, and their simulations """

    class Meta:
        label = 'analyze'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Analyze a random whole-cell model and simulations"
        arguments = [
            (['--config-path'], dict(type=str, default=None, help='Path to configuration file')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        config = rand_wc_model_gen.config.get_config(extra_path=args.config_path)['rand_wc_model_gen']
        kb = wc_kb.io.Reader().run(config['kb']['path']['core'], config['kb']['path']['seq'])
        model = wc_lang.io.Reader().run(config['model']['path'])
        runner = rand_wc_model_gen.analysis.AnalysisRunner(kb, model, config['sim_results']['path'],
                                                           out_path=config['analysis']['path'],
                                                           options=config['analysis'])
        runner.run()


class App(CementApp):
    """ Command line application 

    Attributes:
        results (:obj:`object`): handler results
    """
    class Meta:
        label = 'rand_wc_model_gen'
        base_controller = 'base'
        handlers = [
            BaseController,
            GenerateController,
            SimulateController,
            AnalyzeController,
        ]

    def __init__(self, label=None, **kw):
        """
        Args:
            label (:obj:`str`, optional): label
            **kw (:obj:`dict`, optional)
        """
        self.results = None
        super(App, self).__init__(label=label, **kw)


def main():
    with App() as app:
        app.run()
