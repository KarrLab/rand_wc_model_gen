"""
Takes in a wc_lang model as input, matches RNA/protein IDs with random RNA/protein sequences, molecular weights, and charges, creates transcription, translation, and RNA degradation reactions
for the randomly generated sequences, outputs the updated model into an Excel spreadsheet

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-01-26
:Copyright: 2018, Karr Lab
:License: MIT
"""

from random_wc_model_generator.random_polymer import RandomSeqGen
import wc_lang.io
import wc_lang.core
import xlrd

class Enrich_Polymers (object):

    def __init__ (self, filename): #filename that contains wc_lang model
        #print('enter')
        randomseq = RandomSeqGen()
        core_wc_model = wc_lang.io.Reader().run(filename)
        rnacount = 0
        #Goes through species_types list and finds number of RNAs
        for species_type in core_wc_model.get_species_types():
            if species_type.type == wc_lang.core.SpeciesTypeType.rna:
                rnacount += 1
                
        lists = randomseq.gen_species_types(rnacount) #Uses the RandomSeqGen class to randomly generate RNA/protein sequences (also protein charges and molecular weights) 
        proteins = lists[2]
        seq = proteins[0]
        mw = proteins[1]
        charge = proteins[2]
        rnas = lists[1]

        #Goes through all the RNA IDs in species_types list and matches them with sequence, weight, and charge
        species_types = core_wc_model.get_species_types()
        index = 0
        for species_type in species_types:
            if species_type.type == wc_lang.core.SpeciesTypeType.rna:
                rna = rnas[index]
                species_type.structure = rna
                species_type.molecular_weight = self.rnaweight(rna) #calls on RNA weight calculation function 
                species_type.charge = self.rnacharge(rna) #calls on RNA charge calculation function
                index += 1

        index = 0
        #Goes through all protein IDs in species_types list and matches them with sequence, weight, and charge corresponding to its RNA
        for species_type in species_types:
            if species_type.type == wc_lang.core.SpeciesTypeType.protein:
                species_type.structure = seq[index]
                species_type.molecular_weight = mw[index]
                species_type.charge = charge[index]
                index += 1
            
        self.transcription(core_wc_model, species_types)
        self.translation(core_wc_model, species_types)
        self.rnadegradation(core_wc_model, species_types)
        model_filename = "random_model.xlsx" #outputs updated model to Excel file
        wc_lang.io.Writer().run(model_filename, core_wc_model)


    def rnaweight (self, rna):
        """ Calculates molecular weight for RNA molecule

        Args:
            rna (:obj:`string`): RNA sequence

        Returns:
            :obj:`float`: molecular weight 
        """
        dictionary = {'A':329.2, 'U':306.2, 'C':305.2, 'G':345.2} #Each nucleotide base mapped to its molecular weight (including sugar and phosphate)
        rnamw = 159 #weight of 5' phosphate
        for j in range(len(rna)):
            base = rna[j]
            rnamw += dictionary[base]

        return rnamw

    def rnacharge (self, rna): 
        """ Calculates charge for RNA molecule

        Args:
            rnas (:obj:`string`): RNA sequence

        Returns:
            :obj:`int`: charge of RNA molecule
        """
        return    (-1 * (len(rna) + 1)) #Charge of RNA molecule is the negative of its (nucleotide length + 1)


    def transcription (self, core_model, species_types):
        """
        Creates transcription reactions for the randomly generated RNA sequences

        Args:
            core_model (:obj:`Model`): wc_lang random model

        """
        b = xlrd.open_workbook("Model.xlsx")
        reactions = b.sheet_by_index(3)
        nameList = []
        idList = []
        for i in range(reactions.nrows):
            current_row = reactions.row_values(i)
            if current_row[2] == 'Transcription':
                nameList.append(current_row[1])
                idList.append(current_row[0].replace('-', '_'))
        #write the transcription reactions into the model in memory
        index = 0
        submodels = core_model.get_submodels()
        for submodel in submodels:
            if submodel.name == "Transcription":
                transSubModel = submodel
        for element in species_types:
            if element.type == wc_lang.core.SpeciesTypeType.rna:
                reaction = nameList[index]
                current_id = idList[index]
                transcript = wc_lang.core.Reaction(id=current_id, name = reaction, submodel = transSubModel)
                transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id='ATP'),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=-1 * element.structure.count('A'))
                transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id='GTP'),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=-1 * element.structure.count('G'))
                transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id='UTP'),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=-1 * element.structure.count('U'))
                transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id='CTP'),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=-1 * element.structure.count('C'))
                transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id='H2O'),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=-1)
                transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id='H'),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=1)
                transcript.participants.create(species = wc_lang.core.Species(species_type = element,
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=1)
                transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id='PPI'),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                                              coefficient=element.structure.count('A')+element.structure.count('G')+element.structure.count('U')+element.structure.count('C') )
                index += 1


    def translation (self, core_model, species_types):
        """
        Creates RNA degradation reactions for the randomly generated RNA sequences

        Args:
            core_model (:obj:`Model`): wc_lang random model


        """
        b = xlrd.open_workbook("Model.xlsx")
        reactions = b.sheet_by_index(3)
        nameList = []
        idList = []
        for i in range(reactions.nrows):
            current_row = reactions.row_values(i)
            if current_row[2] == 'Translation':
                nameList.append(current_row[1])
                idList.append(current_row[0].replace('-', '_'))
        #write the translation reactions into the model in memory
        index = 0
        submodels = core_model.get_submodels()
        for submodel in submodels:
            if submodel.name == "Translation":
                transSubModel = submodel
        dictionary ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', 
        'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    
        'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    
        'G':'GLY', 'P':'PRO', 'C':'CYS'}
        gtp = wc_lang.core.Species(species_type = core_model.species_types.get(id='GTP'),
                                                                              compartment = core_model.compartments.get(id='c'))
        h2o = wc_lang.core.Species(species_type = core_model.species_types.get(id='H2O'),
                                                                              compartment = core_model.compartments.get(id='c'))
        gdp = wc_lang.core.Species(species_type = core_model.species_types.get(id='GDP'),
                                                                              compartment = core_model.compartments.get(id='c'))

        pi = wc_lang.core.Species(species_type = core_model.species_types.get(id='PI'),
                                                                              compartment = core_model.compartments.get(id='c'))

        h = wc_lang.core.Species(species_type = core_model.species_types.get(id='H'),
                                                                              compartment = core_model.compartments.get(id='c'))
        
        
        for element in species_types:
            if element.type == wc_lang.core.SpeciesTypeType.protein:
                reaction = nameList[index]
                current_id = idList[index]
                transcript = wc_lang.core.Reaction(id=current_id, name = reaction, submodel = transSubModel)
                for aa in dictionary:
                    count = element.structure.count(aa)
                    if count != 0:
                        transcript.participants.create(species = wc_lang.core.Species(species_type = core_model.species_types.get(id=dictionary[aa]),
                                                                              compartment = core_model.compartments.get(id='c')),
                                                       coefficient = -count)
                length = len(element.structure)
                transcript.participants.create(species = gtp,  coefficient=-1 * (2 * length + 3)  )
                transcript.participants.create(species = h2o, coefficient=-1 * (length + 4))
                transcript.participants.create(species = gdp, coefficient= 2 * length + 3)
                transcript.participants.create(species = pi, coefficient = 2*length + 3)
                transcript.participants.create(species = h , coefficient = 2 * length + 3)
                transcript.participants.create(species = wc_lang.core.Species(species_type = element,
                                                                              compartment = core_model.compartments.get(id='c')), coefficient = 1)
                        
                index += 1

    def rnadegradation (self, core_model, species_types):
        """ Creates RNA degradation reactions for the randomly generated RNA sequences

        Args:
            core_model (:obj:`Model`): wc_lang random model

        """
        b = xlrd.open_workbook("Model.xlsx")
        reactions = b.sheet_by_index(3)
        nameList = []
        idList = []
        for i in range(reactions.nrows):
            current_row = reactions.row_values(i)
            if current_row[2] == 'RnaDegradation':
                nameList.append(current_row[1])
                idList.append(current_row[0].replace('-', '_'))
        #write the rna degradation reactions into the model in memory
        index = 0
        submodels = core_model.get_submodels()
        for submodel in submodels:
            if submodel.name == "RNA degradation":
                transSubModel = submodel


        h2o = wc_lang.core.Species(species_type = core_model.species_types.get(id='H2O'),
                                                                              compartment = core_model.compartments.get(id='c'))
        gmp = wc_lang.core.Species(species_type = core_model.species_types.get(id='GMP'),
                                                                              compartment = core_model.compartments.get(id='c'))

        amp = wc_lang.core.Species(species_type = core_model.species_types.get(id='AMP'),
                                                                              compartment = core_model.compartments.get(id='c'))

        cmp = wc_lang.core.Species(species_type = core_model.species_types.get(id='CMP'),
                                                                              compartment = core_model.compartments.get(id='c'))
        
        ump = wc_lang.core.Species(species_type = core_model.species_types.get(id='UMP'),
                                                                              compartment = core_model.compartments.get(id='c'))

        h = wc_lang.core.Species(species_type = core_model.species_types.get(id='H'),
                                                                              compartment = core_model.compartments.get(id='c'))
        
        
        for element in species_types:
            if element.type == wc_lang.core.SpeciesTypeType.rna:
                reaction = nameList[index]
                current_id = idList[index]
                transcript = wc_lang.core.Reaction(id=current_id, name = reaction, submodel = transSubModel)
                
                transcript.participants.create(species = gmp,  coefficient= element.structure.count('G') )
                transcript.participants.create(species = amp, coefficient= element.structure.count('A')  )
                transcript.participants.create(species = cmp, coefficient= element.structure.count('C')  )
                transcript.participants.create(species = ump, coefficient = element.structure.count('U') )
                transcript.participants.create(species = h2o , coefficient = -1 * ( element.structure.count('A')+element.structure.count('G')+element.structure.count('U')+element.structure.count('C') - 1))
                transcript.participants.create(species = h , coefficient = element.structure.count('A')+element.structure.count('G')+element.structure.count('U')+element.structure.count('C') - 1)
                transcript.participants.create(species = wc_lang.core.Species(species_type = element, compartment = core_model.compartments.get(id='c')), coefficient = -1)
                        
                index += 1

        

