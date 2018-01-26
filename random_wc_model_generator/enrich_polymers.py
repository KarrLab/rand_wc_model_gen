# load model_core
for species_type in model.get_species():
    if species_type.type in [rna]:
        # add data to species_type
        species_type.structure = # from rand code
        species_type.molecular_weight = # from rand code
        species_type.charge = # from rand code

# add other components to the model
# reactions
# Transccription, translation

model_filename = #  random_model.xlsx
wc_lang.io.Writer().run(model_filename, model)

        