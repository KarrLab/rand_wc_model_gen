""" Command line programs for generating random whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-05-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
import rand_wc_model_gen


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


class App(CementApp):
    """ Command line application """
    class Meta:
        label = 'rand_wc_model_gen'
        base_controller = 'base'
        handlers = [
            BaseController,
        ]


def main():
    with App() as app:
        app.run()
