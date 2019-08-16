""" Classes to generate random wc models

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-08-13
:Copyright: 2019, Karr Lab
:License: MIT
"""

import math
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot
import numpy
import os
import pkg_resources
import scipy.constants
import wc_lang
import wc_lang.io
from wc_onto import onto as wc_ontology
from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
from wc_utils.util.units import unit_registry

class RandModelGen(object):
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
        submodel = model.submodels.create(id='submodel_rna')

        # compartment:
        # cytosol

        # volume: 50 al
        c_init_volume  = wc_lang.InitVolume(distribution=wc_ontology['WC:normal_distribution'], mean=50 * 1E-18, std=0.0)
        c_ph = wc_lang.Ph(distribution=wc_ontology['WC:normal_distribution'], mean=7.75, std=0.775)
        c = model.compartments.create(id='c', name='Cytosol', init_volume=c_init_volume, ph=c_ph)
        c.init_density = model.parameters.create(id='density_c', value=1100., 
                                                 units=unit_registry.parse_units('g l^-1'))
        volume_c = model.functions.create(id='volume_c', units=unit_registry.parse_units('l'))

        volume_c.expression, error = wc_lang.FunctionExpression.deserialize(
            f'{c.id} / {c.init_density.id}', 
            self.get_rate_law_context(model))
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
        h2o_structure = wc_lang.ChemicalStructure(value='O', format=wc_lang.ChemicalStructureFormat.SMILES)
        h2o_structure.empirical_formula = OpenBabelUtils.get_formula(h2o_structure.get_structure())
        h2o_structure.molecular_weight = h2o_structure.empirical_formula.get_molecular_weight()
        h2o_structure.charge = h2o_structure.get_structure().GetTotalCharge()
        h2o = model.species_types.create(id='h2o', name='H2O', type=wc_ontology['WC:metabolite'], structure=h2o_structure)
        init_concs['h2o'] = 0.0005 * Avogadro.value * c.init_volume.mean

        h_structure = wc_lang.ChemicalStructure(value='[H+]', format=wc_lang.ChemicalStructureFormat.SMILES)
        h_structure.empirical_formula = OpenBabelUtils.get_formula(h_structure.get_structure())
        h_structure.molecular_weight = h_structure.empirical_formula.get_molecular_weight()
        h_structure.charge = h_structure.get_structure().GetTotalCharge()
        h = model.species_types.create(id='h', name='H', type=wc_ontology['WC:metabolite'], structure=h_structure)
        init_concs['h'] = 0.00005 * Avogadro.value * c.init_volume.mean

        ppi_structure = wc_lang.ChemicalStructure(value='OP(=O)([O-])OP(=O)([O-])[O-]', format=wc_lang.ChemicalStructureFormat.SMILES)
        ppi_structure.empirical_formula = OpenBabelUtils.get_formula(ppi_structure.get_structure())
        ppi_structure.molecular_weight = ppi_structure.empirical_formula.get_molecular_weight()
        ppi_structure.charge = ppi_structure.get_structure().GetTotalCharge()
        ppi = model.species_types.create(id='ppi', name='PPi', type=wc_ontology['WC:metabolite'], structure=ppi_structure)
        init_concs['ppi'] = 0.00005 * Avogadro.value * c.init_volume.mean

        # ntp
        atp_structure = wc_lang.ChemicalStructure(value='C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O)N', format=wc_lang.ChemicalStructureFormat.SMILES)
        atp_structure.empirical_formula = OpenBabelUtils.get_formula(atp_structure.get_structure())
        atp_structure.molecular_weight = atp_structure.empirical_formula.get_molecular_weight()
        atp_structure.charge = atp_structure.get_structure().GetTotalCharge()
        atp = model.species_types.create(id='atp', name='ATP', type=wc_ontology['WC:metabolite'], structure=atp_structure)
        init_concs['atp'] = 0.005 * Avogadro.value * c.init_volume.mean

        gtp_structure = wc_lang.ChemicalStructure(value='C1=NC2=C(N1C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O)N=C(NC2=O)N', format=wc_lang.ChemicalStructureFormat.SMILES)
        gtp_structure.empirical_formula = OpenBabelUtils.get_formula(gtp_structure.get_structure())
        gtp_structure.molecular_weight = gtp_structure.empirical_formula.get_molecular_weight()
        gtp_structure.charge = gtp_structure.get_structure().GetTotalCharge()
        gtp = model.species_types.create(id='gtp', name='GTP', type=wc_ontology['WC:metabolite'], structure=gtp_structure)
        init_concs['gtp'] = 0.005 * Avogadro.value * c.init_volume.mean

        ctp_structure = wc_lang.ChemicalStructure(value='C1=CN(C(=O)N=C1N)C2C(C(C(O2)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O', format=wc_lang.ChemicalStructureFormat.SMILES)
        ctp_structure.empirical_formula = OpenBabelUtils.get_formula(ctp_structure.get_structure())
        ctp_structure.molecular_weight = ctp_structure.empirical_formula.get_molecular_weight()
        ctp_structure.charge = ctp_structure.get_structure().GetTotalCharge()
        ctp = model.species_types.create(id='ctp', name='CTP', type=wc_ontology['WC:metabolite'], structure=ctp_structure)
        init_concs['ctp'] = 0.005 * Avogadro.value * c.init_volume.mean

        utp_structure = wc_lang.ChemicalStructure(value='C1=CN(C(=O)NC1=O)C2C(C(C(O2)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O', format=wc_lang.ChemicalStructureFormat.SMILES)
        utp_structure.empirical_formula = OpenBabelUtils.get_formula(utp_structure.get_structure())
        utp_structure.molecular_weight = utp_structure.empirical_formula.get_molecular_weight()
        utp_structure.charge = utp_structure.get_structure().GetTotalCharge()
        utp = model.species_types.create(id='utp', name='UTP', type=wc_ontology['WC:metabolite'], structure=utp_structure)
        init_concs['utp'] = 0.005 * Avogadro.value * c.init_volume.mean

        # nmp
        amp_structure = wc_lang.ChemicalStructure(value='C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)([O-])[O-])O)O)N', format=wc_lang.ChemicalStructureFormat.SMILES)
        amp_structure.empirical_formula = OpenBabelUtils.get_formula(amp_structure.get_structure())
        amp_structure.molecular_weight = amp_structure.empirical_formula.get_molecular_weight()
        amp_structure.charge = amp_structure.get_structure().GetTotalCharge()
        amp = model.species_types.create(id='amp', name='AMP', type=wc_ontology['WC:metabolite'], structure=amp_structure)
        init_concs['amp'] = 0.00005 * Avogadro.value * c.init_volume.mean

        gmp_structure = wc_lang.ChemicalStructure(value='C1=NC2=C(N1C3C(C(C(O3)COP(=O)([O-])[O-])O)O)N=C(NC2=O)N', format=wc_lang.ChemicalStructureFormat.SMILES)
        gmp_structure.empirical_formula = OpenBabelUtils.get_formula(gmp_structure.get_structure())
        gmp_structure.molecular_weight = gmp_structure.empirical_formula.get_molecular_weight()
        gmp_structure.charge = gmp_structure.get_structure().GetTotalCharge()
        gmp = model.species_types.create(id='gmp', name='GMP', type=wc_ontology['WC:metabolite'], structure=gmp_structure)
        init_concs['gmp'] = 0.00005 * Avogadro.value * c.init_volume.mean

        cmp_structure = wc_lang.ChemicalStructure(value='C1=CN(C(=O)N=C1N)C2C(C(C(O2)COP(=O)([O-])[O-])O)O', format=wc_lang.ChemicalStructureFormat.SMILES)
        cmp_structure.empirical_formula = OpenBabelUtils.get_formula(cmp_structure.get_structure())
        cmp_structure.molecular_weight = cmp_structure.empirical_formula.get_molecular_weight()
        cmp_structure.charge = cmp_structure.get_structure().GetTotalCharge()
        cmp = model.species_types.create(id='cmp', name='CMP', type=wc_ontology['WC:metabolite'], structure=cmp_structure)
        init_concs['cmp'] = 0.00005 * Avogadro.value * c.init_volume.mean

        ump_structure = wc_lang.ChemicalStructure(value='C1=CN(C(=O)NC1=O)C2C(C(C(O2)COP(=O)([O-])[O-])O)O', format=wc_lang.ChemicalStructureFormat.SMILES)
        ump_structure.empirical_formula = OpenBabelUtils.get_formula(ump_structure.get_structure())
        ump_structure.molecular_weight = ump_structure.empirical_formula.get_molecular_weight()
        ump_structure.charge = ump_structure.get_structure().GetTotalCharge()
        ump = model.species_types.create(id='ump', name='UMP', type=wc_ontology['WC:metabolite'], structure=ump_structure)
        init_concs['ump'] = 0.00005 * Avogadro.value * c.init_volume.mean


        # RNA
        # half life = 3 min
        rna_1_str = 'AAUGUGC'
        rna_1_structure = wc_lang.ChemicalStructure(
                            value=rna_1_str,
                            format=wc_lang.ChemicalStructureFormat.BpForms,
                            alphabet=wc_lang.ChemicalStructureAlphabet.rna)
        rna_1_structure.empirical_formula = rna_1_structure.get_structure().get_formula()
        rna_1_structure.molecular_weight = rna_1_structure.empirical_formula.get_molecular_weight()
        rna_1_structure.charge = rna_1_structure.get_structure().get_charge()
        rna_1 = model.species_types.create(id='rna_1',
                                    name='RNA 1',
                                    type=wc_ontology['WC:RNA'],
                                    structure=rna_1_structure)
        half_life_rna_1 = model.parameters.create(id='half_life_rna_1',
                                                  type=None,
                                                  value=180,
                                                  units=unit_registry.parse_units('s'))
        init_concs['rna_1'] = 1

        rna_2_str = 'UCAG'
        rna_2_structure = wc_lang.ChemicalStructure(
                            value=rna_2_str,
                            format=wc_lang.ChemicalStructureFormat.BpForms,
                            alphabet=wc_lang.ChemicalStructureAlphabet.rna)
        rna_2_structure.empirical_formula = rna_2_structure.get_structure().get_formula()
        rna_2_structure.molecular_weight = rna_2_structure.empirical_formula.get_molecular_weight()
        rna_2_structure.charge = rna_2_structure.get_structure().get_charge()
        rna_2 = model.species_types.create(id='rna_2',
                                    name='RNA 2',
                                    type=wc_ontology['WC:RNA'],
                                    structure=rna_2_structure)
        half_life_rna_2 = model.parameters.create(id='half_life_rna_2',
                                                  type=None,
                                                  value=180,
                                                  units=unit_registry.parse_units('s'))
        init_concs['rna_2'] = 1

        rna_3_str = 'ACGUC'
        rna_3_structure = wc_lang.ChemicalStructure(
                            value=rna_3_str,
                            format=wc_lang.ChemicalStructureFormat.BpForms,
                            alphabet=wc_lang.ChemicalStructureAlphabet.rna)
        rna_3_structure.empirical_formula = rna_3_structure.get_structure().get_formula()
        rna_3_structure.molecular_weight = rna_3_structure.empirical_formula.get_molecular_weight()
        rna_3_structure.charge = rna_3_structure.get_structure().get_charge()
        rna_3 = model.species_types.create(id='rna_3',
                                    name='RNA 3',
                                    type=wc_ontology['WC:RNA'],
                                    structure=rna_3_structure)
        half_life_rna_3 = model.parameters.create(id='half_life_rna_3',
                                                  type=None,
                                                  value=180,
                                                  units=unit_registry.parse_units('s'))
        init_concs['rna_3'] = 1

        # enzymes
        rna_pol = model.species_types.create(id='rna_pol', name='RNA polymerase', type=wc_ontology['WC:metabolite'])
        init_concs['rna_pol'] = 10 ** 2
        rna_se = model.species_types.create(id='rna_se', name='RNAse', type=wc_ontology['WC:metabolite'])
        init_concs['rna_se'] = 10 ** 2
        atp_synthase = model.species_types.create(
            id='atp_synthase', 
            name='ATP synthase', 
            type=wc_ontology['WC:protein'])
        init_concs['atp_synthase'] = 10 ** 3
        gtp_synthase = model.species_types.create(
            id='gtp_synthase', 
            name='GTP synthase', 
            type=wc_ontology['WC:protein'])
        init_concs['gtp_synthase'] = 10 ** 3
        ctp_synthase = model.species_types.create(
            id='ctp_synthase', 
            name='CTP synthase', 
            type=wc_ontology['WC:protein'])
        init_concs['ctp_synthase'] = 10 ** 3
        utp_synthase = model.species_types.create(
            id='utp_synthase', 
            name='UTP synthase', 
            type=wc_ontology['WC:protein'])
        init_concs['utp_synthase'] = 10 ** 3

        # species and initial concentrations
        for model_species_type in model.species_types:
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=c)
            model_species.id = model_species.gen_id()
            conc = model.distribution_init_concentrations.create(species=model_species, mean=init_concs[model_species_type.id], units=unit_registry.parse_units('molecule'))
            conc.id = conc.gen_id()

        all_species = {species.id: species for species in model.species}

        # reactions and rate laws


        # rna synthesis (transcription)
        km_atp_trans = model.parameters.create(id='km_atp_trans', value=0.005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        km_gtp_trans = model.parameters.create(id='km_gtp_trans', value=0.005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        km_ctp_trans = model.parameters.create(id='km_ctp_trans', value=0.005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        km_utp_trans = model.parameters.create(id='km_utp_trans', value=0.005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))

        trans_rna_1 = model.reactions.get_or_create(submodel=submodel, id='transcription_' + 'rna_1')
        trans_rna_1.name = 'transcription '+'RNA 1'
        trans_rna_1.participants = []
        # lhs
        trans_rna_1.participants.add(atp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_1_str.count('A')))
        trans_rna_1.participants.add(ctp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_1_str.count('C')))
        trans_rna_1.participants.add(gtp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_1_str.count('G')))
        trans_rna_1.participants.add(utp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_1_str.count('U')))
        trans_rna_1.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        # rhs
        trans_rna_1.participants.add(rna_1.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        trans_rna_1.participants.add(ppi.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=len(rna_1_str)))
        trans_rna_1.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        # rate law
        k_trans_rna_1 = model.parameters.create(
            id='k_trans_rna_1', 
            value=math.log(2)/half_life_rna_1.value * 8, 
            type=wc_ontology['WC:k_cat'], 
            units=unit_registry.parse_units('s^-1 / M'))
        trans_rna_1_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_trans_rna_1'
            ' * (atp[c] / (km_atp_trans * Avogadro * volume_c + atp[c]))'
            ' * (gtp[c] / (km_gtp_trans * Avogadro * volume_c + gtp[c]))'
            ' * (ctp[c] / (km_ctp_trans * Avogadro * volume_c + ctp[c]))'
            ' * (utp[c] / (km_utp_trans * Avogadro * volume_c + utp[c]))'
            ' * rna_pol[c] / (Avogadro * volume_c)',
            self.get_rate_law_context(model))
        trans_rna_1_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=trans_rna_1_rate_law_exp,
                              reaction=trans_rna_1,
                              )
        trans_rna_1_rate_law.id = trans_rna_1_rate_law.gen_id()

        trans_rna_2 = model.reactions.get_or_create(submodel=submodel, id='transcription_' + 'rna_2')
        trans_rna_2.name = 'transcription '+'RNA 2'
        trans_rna_2.participants = []
        # lhs
        trans_rna_2.participants.add(atp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_2_str.count('A')))
        trans_rna_2.participants.add(ctp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_2_str.count('C')))
        trans_rna_2.participants.add(gtp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_2_str.count('G')))
        trans_rna_2.participants.add(utp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_2_str.count('U')))
        trans_rna_2.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        # rhs
        trans_rna_2.participants.add(rna_2.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        trans_rna_2.participants.add(ppi.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=len(rna_2_str)))
        trans_rna_2.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        # rate law
        k_trans_rna_2 = model.parameters.create(
            id='k_trans_rna_2', 
            value=math.log(2)/half_life_rna_2.value * 8, 
            type=wc_ontology['WC:k_cat'], 
            units=unit_registry.parse_units('s^-1 / M'))
        trans_rna_2_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_trans_rna_2'
            ' * (atp[c] / (km_atp_trans * Avogadro * volume_c + atp[c]))'
            ' * (gtp[c] / (km_gtp_trans * Avogadro * volume_c + gtp[c]))'
            ' * (ctp[c] / (km_ctp_trans * Avogadro * volume_c + ctp[c]))'
            ' * (utp[c] / (km_utp_trans * Avogadro * volume_c + utp[c]))'
            ' * rna_pol[c] / (Avogadro * volume_c)',
            self.get_rate_law_context(model))
        trans_rna_2_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=trans_rna_2_rate_law_exp,
                              reaction=trans_rna_2,
                              )
        trans_rna_2_rate_law.id = trans_rna_2_rate_law.gen_id()

        trans_rna_3 = model.reactions.get_or_create(submodel=submodel, id='transcription_' + 'rna_3')
        trans_rna_3.name = 'transcription '+'RNA 3'
        trans_rna_3.participants = []
        # lhs
        trans_rna_3.participants.add(atp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_3_str.count('A')))
        trans_rna_3.participants.add(ctp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_3_str.count('C')))
        trans_rna_3.participants.add(gtp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_3_str.count('G')))
        trans_rna_3.participants.add(utp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-rna_3_str.count('U')))
        trans_rna_3.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        # rhs
        trans_rna_3.participants.add(rna_3.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        trans_rna_3.participants.add(ppi.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=len(rna_3_str)))
        trans_rna_3.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        # rate law
        k_trans_rna_3 = model.parameters.create(
            id='k_trans_rna_3', 
            value=math.log(2)/half_life_rna_3.value * 8, 
            type=wc_ontology['WC:k_cat'], 
            units=unit_registry.parse_units('s^-1 / M'))
        trans_rna_3_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_trans_rna_3'
            ' * (atp[c] / (km_atp_trans * Avogadro * volume_c + atp[c]))'
            ' * (gtp[c] / (km_gtp_trans * Avogadro * volume_c + gtp[c]))'
            ' * (ctp[c] / (km_ctp_trans * Avogadro * volume_c + ctp[c]))'
            ' * (utp[c] / (km_utp_trans * Avogadro * volume_c + utp[c]))'
            ' * rna_pol[c] / (Avogadro * volume_c)',
            self.get_rate_law_context(model))
        trans_rna_3_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=trans_rna_3_rate_law_exp,
                              reaction=trans_rna_3,
                              )
        trans_rna_3_rate_law.id = trans_rna_3_rate_law.gen_id()


        # rna degradation
        deg_rna_1 = model.reactions.get_or_create(submodel=submodel, id='degradation_' + 'rna_1')
        deg_rna_1.name = 'degradation ' + 'RNA 1'
        deg_rna_1.participants = []
        # lhs
        deg_rna_1.participants.add(rna_1.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        deg_rna_1.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-(len(rna_1_str)-1)))
        # rhs
        deg_rna_1.participants.add(amp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_1_str.count('A')))
        deg_rna_1.participants.add(cmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_1_str.count('C')))
        deg_rna_1.participants.add(gmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_1_str.count('G')))
        deg_rna_1.participants.add(ump.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_1_str.count('U')))
        deg_rna_1.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=len(rna_1_str)-1))
        # rate law
        k_deg_rna_1 = model.parameters.create(
            id='k_deg_rna_1', 
            value=math.log(2)/half_life_rna_1.value, 
            type=wc_ontology['WC:k_cat'], 
            units=unit_registry.parse_units('s^-1 / M'))
        km_deg_rna_1 = model.parameters.create(id='km_deg_rna_1', value=1 / Avogadro.value / c.init_volume.mean, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        deg_rna_1_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_deg_rna_1'
            ' * rna_1[c] / (km_deg_rna_1 * Avogadro * volume_c + rna_1[c])'
            ' * rna_se[c] / (Avogadro * volume_c)',
            self.get_rate_law_context(model))
        deg_rna_1_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=deg_rna_1_rate_law_exp,
                              reaction=deg_rna_1,
                              )
        deg_rna_1_rate_law.id = deg_rna_1_rate_law.gen_id()

        deg_rna_2 = model.reactions.get_or_create(submodel=submodel, id='degradation_' + 'rna_2')
        deg_rna_2.name = 'degradation ' + 'RNA 2'
        deg_rna_2.participants = []
        # lhs
        deg_rna_2.participants.add(rna_2.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        deg_rna_2.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-(len(rna_2_str)-1)))
        # rhs
        deg_rna_2.participants.add(amp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_2_str.count('A')))
        deg_rna_2.participants.add(cmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_2_str.count('C')))
        deg_rna_2.participants.add(gmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_2_str.count('G')))
        deg_rna_2.participants.add(ump.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_2_str.count('U')))
        deg_rna_2.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=len(rna_2_str)-1))
        # rate law
        k_deg_rna_2 = model.parameters.create(
            id='k_deg_rna_2', 
            value=math.log(2)/half_life_rna_2.value, 
            type=wc_ontology['WC:k_cat'], 
            units=unit_registry.parse_units('s^-1 / M'))
        km_deg_rna_2 = model.parameters.create(id='km_deg_rna_2', value=1 / Avogadro.value / c.init_volume.mean, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        deg_rna_2_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_deg_rna_2'
            ' * rna_2[c] / (km_deg_rna_2 * Avogadro * volume_c + rna_2[c])'
            ' * rna_se[c] / (Avogadro * volume_c)',
            self.get_rate_law_context(model))
        deg_rna_2_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=deg_rna_2_rate_law_exp,
                              reaction=deg_rna_2,
                              )
        deg_rna_2_rate_law.id = deg_rna_2_rate_law.gen_id()

        deg_rna_3 = model.reactions.get_or_create(submodel=submodel, id='degradation_' + 'rna_3')
        deg_rna_3.name = 'degradation ' + 'RNA 3'
        deg_rna_3.participants = []
        # lhs
        deg_rna_3.participants.add(rna_3.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        deg_rna_3.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-(len(rna_3_str)-1)))
        # rhs
        deg_rna_3.participants.add(amp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_3_str.count('A')))
        deg_rna_3.participants.add(cmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_3_str.count('C')))
        deg_rna_3.participants.add(gmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_3_str.count('G')))
        deg_rna_3.participants.add(ump.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=rna_3_str.count('U')))
        deg_rna_3.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=len(rna_3_str)-1))
        # rate law
        k_deg_rna_3 = model.parameters.create(
            id='k_deg_rna_3', 
            value=math.log(2)/half_life_rna_3.value, 
            type=wc_ontology['WC:k_cat'], 
            units=unit_registry.parse_units('s^-1 / M'))
        km_deg_rna_3 = model.parameters.create(id='km_deg_rna_3', value=1 / Avogadro.value / c.init_volume.mean, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        deg_rna_3_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_deg_rna_3'
            ' * rna_3[c] / (km_deg_rna_3 * Avogadro * volume_c + rna_3[c])'
            ' * rna_se[c] / (Avogadro * volume_c)',
            self.get_rate_law_context(model))
        deg_rna_3_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=deg_rna_3_rate_law_exp,
                              reaction=deg_rna_3,
                              )
        deg_rna_3_rate_law.id = deg_rna_3_rate_law.gen_id()


        # ntp synthesis from nmp
        km_syn_ntp_ppi = model.parameters.create(id='km_syn_ntp_ppi', value=0.00005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))

        # atp
        syn_atp = model.reactions.get_or_create(submodel=submodel, id='syn_atp')
        syn_atp.name = 'synthesis ' + 'ATP'
        syn_atp.participants = []
        # lhs
        syn_atp.participants.add(amp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_atp.participants.add(ppi.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_atp.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        # rhs
        syn_atp.participants.add(atp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        syn_atp.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        # rate law
        k_syn_atp = model.parameters.create(id='k_syn_atp', value=math.log(2)/half_life_rna_3.value * 2 * 4, type=wc_ontology['WC:k_cat'], units=unit_registry.parse_units('s^-1'))
        km_syn_atp_amp = model.parameters.create(id='km_syn_atp_amp', value=0.00005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        syn_atp_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_syn_atp * (amp[c] / (km_syn_atp_amp * Avogadro * volume_c + amp[c])) * (ppi[c] / (km_syn_ntp_ppi * Avogadro * volume_c + ppi[c]))',
            self.get_rate_law_context(model))
        syn_atp_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=syn_atp_rate_law_exp,
                              reaction=syn_atp,
                              )
        syn_atp_rate_law.id = syn_atp_rate_law.gen_id()

        # gtp
        syn_gtp = model.reactions.get_or_create(submodel=submodel, id='syn_gtp')
        syn_gtp.name = 'synthesis ' + 'GTP'
        syn_gtp.participants = []
        # lhs
        syn_gtp.participants.add(gmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_gtp.participants.add(ppi.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_gtp.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        # rhs
        syn_gtp.participants.add(gtp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        syn_gtp.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        # rate law
        k_syn_gtp = model.parameters.create(id='k_syn_gtp', value=math.log(2)/half_life_rna_3.value * 2 * 4, type=wc_ontology['WC:k_cat'], units=unit_registry.parse_units('s^-1'))
        km_syn_gtp_gmp = model.parameters.create(id='km_syn_gtp_gmp', value=0.00005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        syn_gtp_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_syn_gtp * (gmp[c] / (km_syn_gtp_gmp * Avogadro * volume_c + gmp[c])) * (ppi[c] / (km_syn_ntp_ppi * Avogadro * volume_c + ppi[c]))',
            self.get_rate_law_context(model))
        syn_gtp_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=syn_gtp_rate_law_exp,
                              reaction=syn_gtp,
                              )
        syn_gtp_rate_law.id = syn_gtp_rate_law.gen_id()

        # ctp
        syn_ctp = model.reactions.get_or_create(submodel=submodel, id='syn_ctp')
        syn_ctp.name = 'synthesis ' + 'CTP'
        syn_ctp.participants = []
        # lhs
        syn_ctp.participants.add(cmp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_ctp.participants.add(ppi.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_ctp.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        # rhs
        syn_ctp.participants.add(ctp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        syn_ctp.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        # rate law
        k_syn_ctp = model.parameters.create(id='k_syn_ctp', value=math.log(2)/half_life_rna_3.value * 2 * 4, type=wc_ontology['WC:k_cat'], units=unit_registry.parse_units('s^-1'))
        km_syn_ctp_cmp = model.parameters.create(id='km_syn_ctp_cmp', value=0.00005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        syn_ctp_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_syn_ctp * (cmp[c] / (km_syn_ctp_cmp * Avogadro * volume_c + cmp[c])) * (ppi[c] / (km_syn_ntp_ppi * Avogadro * volume_c + ppi[c]))',
            self.get_rate_law_context(model))
        syn_ctp_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=syn_ctp_rate_law_exp,
                              reaction=syn_ctp,
                              )
        syn_ctp_rate_law.id = syn_ctp_rate_law.gen_id()

        # utp
        syn_utp = model.reactions.get_or_create(submodel=submodel, id='syn_utp')
        syn_utp.name = 'synthesis ' + 'UTP'
        syn_utp.participants = []
        # lhs
        syn_utp.participants.add(ump.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_utp.participants.add(ppi.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        syn_utp.participants.add(h.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=-1))
        # rhs
        syn_utp.participants.add(utp.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        syn_utp.participants.add(h2o.species.get_one(compartment=c).species_coefficients.get_or_create(coefficient=1))
        # rate law
        k_syn_utp = model.parameters.create(id='k_syn_utp', value=math.log(2)/half_life_rna_3.value * 2 * 4, type=wc_ontology['WC:k_cat'], units=unit_registry.parse_units('s^-1'))
        km_syn_utp_ump = model.parameters.create(id='km_syn_utp_ump', value=0.00005, type=wc_ontology['WC:K_m'], units=unit_registry.parse_units('M'))
        syn_utp_rate_law_exp, errors = wc_lang.RateLawExpression.deserialize(
            'k_syn_utp * (ump[c] / (km_syn_utp_ump * Avogadro * volume_c + ump[c])) * (ppi[c] / (km_syn_ntp_ppi * Avogadro * volume_c + ppi[c]))',
            self.get_rate_law_context(model))
        syn_utp_rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                              type=None,
                              expression=syn_utp_rate_law_exp,
                              reaction=syn_utp,
                              )
        syn_utp_rate_law.id = syn_utp_rate_law.gen_id()


        return model

    @classmethod
    def get_rate_law_context(cls, model):
        return {
            wc_lang.Compartment: cls.get_rate_law_compartment_context(model),
            wc_lang.Species: cls.get_rate_law_species_context(model),
            wc_lang.Parameter: cls.get_rate_law_parameter_context(model),            
            wc_lang.Function: cls.get_rate_law_function_context(model),
        }

    @classmethod
    def get_rate_law_compartment_context(cls, model):
        return {compartment.id: compartment for compartment in model.compartments}

    @classmethod
    def get_rate_law_species_context(cls, model):
        return {species.id: species for species in model.species}

    @classmethod
    def get_rate_law_parameter_context(cls, model):
        return {parameter.id: parameter for parameter in model.parameters}

    @classmethod
    def get_rate_law_function_context(cls, model):
        return {function.id: function for function in model.functions}


if __name__ == '__main__':

    model_filename = pkg_resources.resource_filename('rand_wc_model_gen', os.path.join('model_gen', 'model.xlsx'))
    results_parent_dirname = 'results'
    checkpoint_period = 100.
    end_time = 100.

    # generate model
    model = RandModelGen(options={'id':'test_rand', 'name':'test random model', 'version':'0.0'}).run()

    # write model
    wc_lang.io.Writer().run(model_filename, model, data_repo_metadata=False)

    # model = wc_lang.io.Reader().run(model_filename)[wc_lang.Model][0]
    
    # simulate model
    sim = Simulation(model)
    _, results_dirname = sim.run(end_time=end_time,
                                 results_dir=results_parent_dirname,
                                 checkpoint_period=checkpoint_period)
    results = RunResults(results_dirname)
