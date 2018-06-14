""" Tests of generation of chromosomes and genes

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen.kb_gen import metabolites
import csv
import numpy
import os
import tempfile
import unittest
import wc_kb
import wc_utils.util.chem


class MetabolitesGeneratorTestCase(unittest.TestCase):
    def setUp(self):
        fid, self.data_path = tempfile.mkstemp()
        os.close(fid)

        with open(self.data_path, 'w') as file:
            writer = csv.DictWriter(file, fieldnames=('Id', 'Name', 'Structure (InChI)', 'Charge', 'Intracellular concentration (M)'))
            writer.writeheader()

            mets = [
                {
                    'Id': 'h2o',
                    'Name': 'Water',
                    'Structure (InChI)': 'InChI=1S/H2O/h1H2',
                    'Intracellular concentration (M)': 55.,
                },
                {
                    'Id': 'h',
                    'Name': 'hydrogen(1+)',
                    'Structure (InChI)': 'InChI=1S/p+1',
                    'Intracellular concentration (M)': 1e-6,
                },
            ]
            for met in mets:
                writer.writerow(met)

    def tearDown(self):
        os.remove(self.data_path)

    def test_run(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        gen = metabolites.MetabolitesGenerator(kb, options={
            'data_path': self.data_path,
        })
        gen.run()

        h2o = cell.species_types.get_one(__type=wc_kb.MetaboliteSpeciesType, id='h2o')
        self.assertEqual(h2o.get_empirical_formula(), wc_utils.util.chem.EmpiricalFormula('H2O'))
        self.assertEqual(h2o.get_charge(), 0)
        self.assertEqual(h2o.concentration, 55.)

        h = cell.species_types.get_one(__type=wc_kb.MetaboliteSpeciesType, id='h')
        self.assertEqual(h.get_empirical_formula(), wc_utils.util.chem.EmpiricalFormula('H'))
        self.assertEqual(h.get_charge(), 1)
        self.assertEqual(h.concentration, 1e-6)
