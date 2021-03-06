""" Command line programs for generating random whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

import cement
import os
import rand_wc_model_gen
import rand_wc_model_gen.analysis
import rand_wc_model_gen.config
import wc_kb
import wc_kb.io
import wc_kb_gen.random
import wc_lang
import wc_lang.io
import wc_model_gen.prokaryote
import wc_sim.run_results
import wc_sim.simulation


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "Random whole-cell model generator"
        help = "Random whole-cell model generator"
        arguments = [
            (['-v', '--version'], dict(action='version', version=rand_wc_model_gen.__version__)),
        ]

    @cement.ex(hide=True)
    def _default(self):
        self._parser.print_help()


class GenerateController(cement.Controller):
    """ Generate a random whole-cell knowledge base and a random whole-cell model """

    class Meta:
        label = 'generate'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Generate a random whole-cell knowledge base and a random whole-cell model"
        help = "Generate a random whole-cell knowledge base and a random whole-cell model"
        arguments = [
            (['--config-path'], dict(type=str,
                                     default=None,
                                     help='Path to configuration file')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        config = rand_wc_model_gen.config.get_config(extra_path=args.config_path)['rand_wc_model_gen']
        kb = wc_kb_gen.random.RandomKbGenerator(options=config['kb_gen']).run()
        model = wc_model_gen.prokaryote.ProkaryoteModelGenerator(kb, options=config['model_gen']).run()

        if not os.path.isdir(os.path.dirname(config['kb']['path']['core'])):
            os.makedirs(os.path.dirname(config['kb']['path']['core']))
        if not os.path.isdir(os.path.dirname(config['kb']['path']['seq'])):
            os.makedirs(os.path.dirname(config['kb']['path']['seq']))
        wc_kb.io.Writer().run(config['kb']['path']['core'], kb,
                              seq_path=config['kb']['path']['seq'],
                              data_repo_metadata=config['kb_gen']['data_repo_metadata'])

        if not os.path.isdir(os.path.dirname(config['model']['path'])):
            os.makedirs(os.path.dirname(config['model']['path']))
        wc_lang.io.Writer().run(config['model']['path'], model,
                                data_repo_metadata=config['model_gen']['data_repo_metadata'])


class SimulateController(cement.Controller):
    """ Simulate a random whole-cell model """

    class Meta:
        label = 'simulate'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Simulate a random whole-cell model"
        help = "Simulate a random whole-cell model"
        arguments = [
            (['--config-path'], dict(type=str, default=None, help='Path to configuration file')),
            (['--seed'], dict(type=int, default=None, help='Random number generator seed')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        config = rand_wc_model_gen.config.get_config(extra_path=args.config_path)['rand_wc_model_gen']

        model = wc_lang.io.Reader().run(config['model']['path'])[wc_lang.Model][0]

        seed = args.seed
        if seed is None:
            seed = config['sim']['seed']
        simulation = wc_sim.simulation.Simulation(model)
        num_events, sim_results_path = simulation.run(end_time=config['sim']['end_time'],
                                                      time_step=config['sim']['time_step'],
                                                      seed=seed,
                                                      checkpoint_period=config['sim_results']['checkpoint_period'],
                                                      results_dir=config['sim_results']['path'])
        self.app.results = {
            'num_events': num_events,
            'sim_results_path': sim_results_path,
        }


class AnalyzeController(cement.Controller):
    """ Analyze a random whole-cell knowledge base, model, and their simulations """

    class Meta:
        label = 'analyze'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Analyze a random whole-cell model and simulations"
        help = "Analyze a random whole-cell model and simulations"
        arguments = [
            (['--config-path'], dict(type=str, default=None, help='Path to configuration file')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        config = rand_wc_model_gen.config.get_config(extra_path=args.config_path)['rand_wc_model_gen']
        kb = wc_kb.io.Reader().run(config['kb']['path']['core'], config['kb']['path']['seq'])[wc_kb.KnowledgeBase][0]
        model = wc_lang.io.Reader().run(config['model']['path'])[wc_lang.Model][0]
        runner = rand_wc_model_gen.analysis.AnalysisRunner(kb, model, config['sim_results']['path'],
                                                           out_path=config['analysis']['path'],
                                                           options=config['analysis'])
        runner.run()


class App(cement.App):
    """ Command line application 

    Attributes:
        results (:obj:`object`): handler results
    """
    class Meta:
        label = 'rand-wc-model-gen'
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
