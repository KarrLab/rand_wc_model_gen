""" Command line programs for generating random whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
import random_wc_model_generator


class BaseController(CementBaseController):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "Random whole-cell model generator"
        arguments = [
            (['-v', '--version'], dict(action='version', version=random_wc_model_generator.__version__)),
        ]

    @expose(hide=True)
    def default(self):
        self.app.args.print_help()


class App(CementApp):
    """ Command line application """
    class Meta:
        label = 'random_wc_model_generator'
        base_controller = 'base'
        handlers = [
            BaseController,
        ]


def main():
    with App() as app:
        app.run()
