""" Generate WC models that can be used to test WC tools

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-01-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

import argparse
import yaml
import wc_lang

'''
Design:

Steps:
    Parse command-line arguments
    Read configuration file(s)
    Create some of the microorganism's structure
        Genome, including genes, promoters, and perhaps operons
        Reaction network properties
        Other physical properties derived from the inputs
    Create a WC model's components as wc_lang instances:
        Parameters based on inputs, such as the cell cycle duration, the fraction of cell mass which is not water, and RNA's half life
        Compartments, including at least extracellular space and the cytosol
        Submodels, including metabolism, transcription, translation, protein and RNA degradation, and DNA replication
        Essential species types, such as ATP, ADP, DNA and RNA nucleotides, tRNAs, RNA polymerase, DNA polymerase, and ribosomal RNA, etc.
        Specie types, including metabolites, DNA, RNAs, and proteins
        Initial specie type concentrations in compartments
        Reactions, with rate-laws
    Save generated model to file(s)
    Summarize generated model on standard out

External steps:
    Create inputs
    Run model generator
    Test generated models
    Use generated models
'''


class CreateWcLangModel(object):
    """ Generate a wc_lang representation of a synthetic whole-cell model.

    """

    def __init__(self):
        pass

    def create_parameters(self):
        pass

    def create_compartments(self):
        pass

    def create_submodels(self):
        pass

    def create_essential_specie_types(self):
        pass

    def create_specie_types(self):
        pass

    def create_concentrations(self):
        pass

    def create_reactions(self):
        pass

    def run(self):
        self.create_parameters()
        self.create_compartments()
        self.create_submodels()
        self.create_essential_specie_types()
        self.create_specie_types()
        self.create_concentrations()
        self.create_reactions()


class GenerateModel(object):
    """ Generate a synthetic whole-cell model

    Attributes:
        args (:obj:`Namespace`): command line arguments
        create_wc_lang_model (:obj:`CreateWcLangModel`): object that creates the wc_lang model
        config (:obj:`dict`): configuration from the config file
    """

    def __init__(self, args):
        self.args = args

    @staticmethod
    def parse_command_line(test_args=None):
        """ Parse the command line

        Args:
            test_args (:obj:`list`): arguments for testing parse_command_line

        Returns:
            :obj:`Namespace`: parsed command-line arguments
        """
        parser = argparse.ArgumentParser(description="Generate WC models for testing WC tools")
        parser.add_argument('config_file', type=argparse.FileType('r'),
            help="Configuration file, in YAML")
        parser.add_argument('generated_model', type=argparse.FileType('w'),
            help="Output file for wc_lang model; if model_format is in {'tsv', 'csv'}, "
                "root for set of delimited files")
        parser.add_argument('--model_format', default='xlsx', choices=['xlsx', 'tsv', 'csv'],
            help="Format for wc_lang model file; default: %(default)s")
        if test_args is None:
            return parser.parse_args()  # pragma: no cover
        else:
            return parser.parse_args(test_args)

    def read_config_file(self):
        self.config = yaml.load(self.args.config_file)

    def create_cellular_structure(self):
        pass

    def create_wc_lang_model(self):
        self.create_wc_lang_model = CreateWcLangModel()
        self.create_wc_lang_model.run()

    def save_generated_model(self, wc_lang_model):
        pass

    def summarize_generated_model(self, wc_lang_model):
        pass

    def run(self):
        """ Generate and save a model

        Returns:
            :obj:`GenerateModel`: the GenerateModel object
        """
        self.read_config_file()
        self.create_cellular_structure()
        wc_lang_model = self.create_wc_lang_model()
        self.save_generated_model(wc_lang_model)
        self.summarize_generated_model(wc_lang_model)
        return self

# todo: move to __main__.py
def main(test_args=None):
    if test_args is None:
        args = GenerateModel.parse_command_line()   # pragma: no cover
    else:
        args = GenerateModel.parse_command_line(test_args)
    return GenerateModel(args).run()

if __name__ == '__main__':
    main()  # pragma: no cover
