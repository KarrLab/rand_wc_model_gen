import os

from wc_lang.core import Compartment, Reaction, SpeciesType, Species
from wc_lang.io import Reader, Writer
from obj_model.core import Validator

model_core = os.path.join(os.path.dirname(__file__), '../rand_wc_model_gen/data/fixtures', 'model_core.xlsx')
model = Reader().run(model_core)

AspAsp = model.get_component('species_type', 'AspAsp')
# AspAsp.pprint(max_depth=0)

ASP_c = list(filter(lambda s: s.id() == 'ASP[c]', model.get_species()))
if len(ASP_c) == 1:
    # ASP_c[0].pprint(max_depth=1)
    pass

species_type_list = [SpeciesType(id="species_{}".format(i), model=model) for i in range(3)]
comp = Compartment(id='comp')
species_list = [Species(species_type=st, compartment=comp) for st in species_type_list]
rxn_1 = Reaction(id='test1')
rxn_2 = Reaction(id='test2')
rxn_1.participants.add(species_list[0].species_coefficients.get_or_create(coefficient=1))
rxn_2.participants.add(species_list[0].species_coefficients.get_or_create(coefficient=1))

error = Validator().run(rxn_1.participants + rxn_2.participants)
print(error)
