from rand_wc_model_gen.kb_gen.genome import GenomeGenerator
import wc_kb

kb = wc_kb.KnowledgeBase()
gen = GenomeGenerator(kb, {})
print(gen.options)
gen.clean_and_validate_options()
gen.gen_components()
