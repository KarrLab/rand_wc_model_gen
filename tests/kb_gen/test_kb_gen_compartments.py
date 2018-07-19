import wc_kb
import wc_kb_gen
import unittest
from rand_wc_model_gen.kb_gen import compartments


class CompartmentsGeneratorTestCase(unittest.TestCase):
    def test_run(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        gen = compartments.CompartmentsGenerator(kb, options={'compartments': [
                                                 ('c', 'cytosol'), ('e', 'extracellular space'), ('t', 'test')]})

        gen.run()
        cytosol = cell.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

        extra = cell.compartments.get_one(id='e')
        self.assertEqual(extra.name, 'extracellular space')

        test = cell.compartments.get_one(id='t')
        self.assertEqual(test.name, 'test')

    def test_default(self):
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()
        gen = compartments.CompartmentsGenerator(kb)

        gen.run()
        cytosol = cell.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

        extra = cell.compartments.get_one(id='e')
        self.assertEqual(extra.name, 'extracellular space')

    def test_validation(self):
        kb = wc_kb.KnowledgeBase()

        with self.assertRaisesRegex(AssertionError, "Could not find cytosol"):
            gen = compartments.CompartmentsGenerator(kb, options={'compartments': [
                ('d', 'cytosol'), ('e', 'extracellular space'), ('t', 'test')]})

        with self.assertRaisesRegex(AssertionError, "Could not find extracellular space"):
            gen = compartments.CompartmentsGenerator(kb, options={'compartments': [
                ('c', 'cytosol'), ('d', 'extracellular space'), ('t', 'test')]})

        with self.assertRaisesRegex(AssertionError, "Must define at least 2 compartments"):
            gen = compartments.CompartmentsGenerator(kb, options={'compartments': [
                ('c', 'cytosol')]})
