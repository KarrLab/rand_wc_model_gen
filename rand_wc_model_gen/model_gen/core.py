""" Classes to generate random wc models

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-08-13
:Copyright: 2019, Karr Lab
:License: MIT
"""

import Bio.Alphabet
import Bio.Seq
from matplotlib import pyplot
import numpy
import scipy.constants
import wc_lang
import wc_lang.io
from wc_onto import onto as wc_ontology
from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults
from wc_utils.util.units import unit_registry

class RandomModelGenerator(object):
    """ Generator for random wc models

    Attributes:
        options (:obj:`dict`, optional): dictionary of options

    """

    def __init__(self, options=None):
        """
        Args:
            options (:obj:`dict`, optional): dictionary of options
        """

        self.options = options or {}
        self.clean_and_validate_options()

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        id = options.get('id', None)
        assert(isinstance(id, str) or id is None)
        options['id'] = id

        name = options.get('name', None)
        assert(isinstance(name, str) or name is None)
        options['name'] = name

        version = options.get('version', None)
        assert(isinstance(version, str) or version is None)
        options['version'] = version

    def run(self):
        """ Generate a :obj:`wc_lang` model

        Returns:
            :obj:`wc_lang.Model`: model
        """
        model = wc_lang.Model()
        model.id = self.options.get('id')
        model.name = self.options.get('name')
        model.version = self.options.get('version')

        # environment
        # model.environment = wc_lang.Environment(id='env', temp=37.0, temp_units=unit_registry.parse_units('degC'))

        # minimal rna model
        submodel = wc_lang.Submodel(id='submodel_rna', model=model)

        # compartment:
        # cytosol

        # volume: 50 al
        c_init_volume  = wc_lang.InitVolume(distribution=wc_ontology['WC:normal_distribution'], mean=50 * 1E-18, std=0)
        c_ph = wc_lang.Ph(distribution=wc_ontology['WC:normal_distribution'], mean=7.75, std=0.775)
        c = wc_lang.Compartment(id='c', name='Cytosol', model=model, init_volume=c_init_volume, ph=c_ph)
        c.init_density = model.parameters.create(id='density_c', value=1100., units=unit_registry.parse_units('g l^-1'))
        volume_c = model.functions.create(id='volume_c', units=unit_registry.parse_units('l'))

        volume_c.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
            wc_lang.Compartment: {c.id: c},
            wc_lang.Parameter: {c.init_density.id: c.init_density},
            })
        assert error is None, str(error)

        # parameters:
        Avogadro = model.parameters.create(id='Avogadro',
                                           type=None,
                                           value=scipy.constants.Avogadro,
                                           units=unit_registry.parse_units('molecule mol^-1'))
        cellCycleLength = model.parameters.create(id='cellCycleLength',
                                                  type=None,
                                                  value=21600,
                                                  units=unit_registry.parse_units('s'))
        # species types

        init_concs = {}

        # other
        h2o = wc_lang.SpeciesType(id='h2o', name='H2O', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['h2o'] = 0.0005 * Avogadro.value * c.init_volume.mean
        h = wc_lang.SpeciesType(id='h', name='H', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['h'] = 0.00005 * Avogadro.value * c.init_volume.mean
        ppi = wc_lang.SpeciesType(id='ppi', name='PPi', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['ppi'] = 0.00005 * Avogadro.value * c.init_volume.mean

        # ntp
        atp = wc_lang.SpeciesType(id='atp', name='ATP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['atp'] = 0.005 * Avogadro.value * c.init_volume.mean
        gtp = wc_lang.SpeciesType(id='gtp', name='GTP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['gtp'] = 0.005 * Avogadro.value * c.init_volume.mean
        ctp = wc_lang.SpeciesType(id='ctp', name='CTP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['ctp'] = 0.005 * Avogadro.value * c.init_volume.mean
        utp = wc_lang.SpeciesType(id='utp', name='UTP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['utp'] = 0.005 * Avogadro.value * c.init_volume.mean

        # nmp
        amp = wc_lang.SpeciesType(id='amp', name='AMP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['amp'] = 0.00005 * Avogadro.value * c.init_volume.mean
        gmp = wc_lang.SpeciesType(id='gmp', name='GMP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['gmp'] = 0.00005 * Avogadro.value * c.init_volume.mean
        cmp = wc_lang.SpeciesType(id='cmp', name='CMP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['cmp'] = 0.00005 * Avogadro.value * c.init_volume.mean
        ump = wc_lang.SpeciesType(id='ump', name='UMP', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['ump'] = 0.00005 * Avogadro.value * c.init_volume.mean


        # RNA
        # half life = 3 min
        rna_1 = wc_lang.SpeciesType(id='rna_1',
                                    name='RNA 1',
                                    type=wc_ontology['WC:RNA'],
                                    structure=wc_lang.ChemicalStructure(
                                        value='AAUG',
                                        format=wc_lang.ChemicalStructureFormat(1),
                                        alphabet=wc_lang.ChemicalStructureAlphabet.rna),
                                    model=model)
        half_life_rna_1 = model.parameters.create(id='half_life_rna_1',
                                                  type=None,
                                                  value=180,
                                                  units=unit_registry.parse_units('s'))
        init_concs['rna_1'] = 1

        rna_2 = wc_lang.SpeciesType(id='rna_2',
                                    name='RNA 2',
                                    type=wc_ontology['WC:RNA'],
                                    structure=wc_lang.ChemicalStructure(
                                        value='UCAG',
                                        format=wc_lang.ChemicalStructureFormat(1),
                                        alphabet=wc_lang.ChemicalStructureAlphabet.rna),
                                    model=model)
        half_life_rna_2 = model.parameters.create(id='half_life_rna_2',
                                                  type=None,
                                                  value=180,
                                                  units=unit_registry.parse_units('s'))
        init_concs['rna_2'] = 1

        rna_3 = wc_lang.SpeciesType(id='rna_3',
                                    name='RNA 3',
                                    type=wc_ontology['WC:RNA'],
                                    structure=wc_lang.ChemicalStructure(
                                        value='ACGU',
                                        format=wc_lang.ChemicalStructureFormat(1),
                                        alphabet=wc_lang.ChemicalStructureAlphabet.rna),
                                    model=model)
        half_life_rna_3 = model.parameters.create(id='half_life_rna_3',
                                                  type=None,
                                                  value=180,
                                                  units=unit_registry.parse_units('s'))
        init_concs['rna_3'] = 1

        # enzymes
        rna_pol = wc_lang.SpeciesType(id='rna_pol', name='RNA polymerase', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['rna_pol'] = 10 ** 2
        rna_se = wc_lang.SpeciesType(id='rna_se', name='RNAse', type=wc_ontology['WC:metabolite'], model=model)
        init_concs['rna_se'] = 10 ** 2
        atp_synthase = wc_lang.SpeciesType(id='atp_synthase', name='ATP synthase', type=wc_ontology['WC:protein'], model=model)
        init_concs['atp_synthase'] = 10 ** 3
        gtp_synthase = wc_lang.SpeciesType(id='gtp_synthase', name='GTP synthase', type=wc_ontology['WC:protein'], model=model)
        init_concs['gtp_synthase'] = 10 ** 3
        ctp_synthase = wc_lang.SpeciesType(id='ctp_synthase', name='CTP synthase', type=wc_ontology['WC:protein'], model=model)
        init_concs['ctp_synthase'] = 10 ** 3
        utp_synthase = wc_lang.SpeciesType(id='utp_synthase', name='UTP synthase', type=wc_ontology['WC:protein'], model=model)
        init_concs['utp_synthase'] = 10 ** 3

        # species and initial concentrations
        for model_species_type in model.species_types:
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=c)
            model_species.id = model_species.gen_id()
            conc = model.distribution_init_concentrations.create(species=model_species, mean=init_concs[model_species_type.id], units=unit_registry.parse_units('molecule'))
            conc.id = conc.gen_id()

        # reactions and rate laws


        return model


if __name__ == '__main__':

    model_filename = 'model.xlsx'
    results_parent_dirname = 'results'
    checkpoint_period = 100.
    end_time = 100.

    # generate model
    model = RandomModelGenerator(options={'id':'test_rand', 'name':'test random model', 'version':'0.0'}).run()

    # # write model
    wc_lang.io.Writer().run(model_filename, model, data_repo_metadata=False)

    # simulate model
    # sim = Simulation(model)
    # _, results_dirname = sim.run(end_time=end_time,
    #                              results_dir=results_parent_dirname,
    #                              checkpoint_period=checkpoint_period)
    # results = RunResults(results_dirname)
