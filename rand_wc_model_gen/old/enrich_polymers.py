"""
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-01-26
:Copyright: 2018, Karr Lab
:License: MIT
"""

from rand_wc_model_gen.random_polymer import RandomSeqGen
import wc_lang.io
import wc_lang.core
import os


# TODO(Arthur): catch exceptions
class Enrich_Polymers(object):
    """
    Take in a `wc_lang` model as input, matches RNA/protein IDs with random RNA/protein sequences,
    molecular weights, and charges, creates transcription, translation, and RNA degradation reactions
    for the randomly generated sequences, outputs an updated `wc_lang` model in an Excel spreadsheet
    """
    # Molecular weight of each RNA nucleotide base (including sugar and phosphate)
    BASE_WEIGHTS = {'A':329.2, 'U':306.2, 'C':305.2, 'G':345.2}
    RNAMW = 159 # weight of 5' phosphate
    CORE_POLYMER_RXN_INFO = os.path.join(os.path.dirname(__file__), 'data/fixtures/core_polymer_rxn_info.csv')

    # Map from one-letter to three-letter amino acid abbreviations
    ABBREVIATIONS ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN',
    'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',
    'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',
    'G':'GLY', 'P':'PRO', 'C':'CYS'}

    def __init__(self, filename, model_filename):
        """
        Args:
            filename (:obj:`str`): filename that contains wc_lang model
            model_filename (:obj: 'str'): filename with the updated output model
        """
        randomseq = RandomSeqGen()
        core_wc_model = wc_lang.io.Reader().run(filename)
        rnacount = 0
        # Go through species_types list and finds number of RNAs
        for species_type in core_wc_model.get_species_types():
            if species_type.type == wc_lang.core.SpeciesTypeType.rna:
                rnacount += 1

        # Use the RandomSeqGen class to randomly generate RNA/protein sequences (also protein charges and molecular weights)
        (genome, genes, rnas, proteins) = randomseq.gen_polymers(rnacount)
        (seq, mw, charge) = proteins

        # Go through all the RNA IDs in species_types list and matches them with sequence, weight, and charge
        species_types = core_wc_model.get_species_types()
        index = 0
        for species_type in species_types:
            if species_type.type == wc_lang.core.SpeciesTypeType.rna:
                rna = rnas[index]
                species_type.structure = rna
                species_type.molecular_weight = self.rna_weight(rna)
                species_type.charge = self.rna_charge(rna)
                index += 1

        index = 0
        # Go through all protein IDs in species_types list and matches them with sequence, weight, and charge corresponding to its RNA
        for species_type in species_types:
            if species_type.type == wc_lang.core.SpeciesTypeType.protein:
                species_type.structure = seq[index]
                species_type.molecular_weight = mw[index]
                species_type.charge = charge[index]
                index += 1

        self.transcription(core_wc_model, species_types)
        self.translation(core_wc_model, species_types)
        self.rna_degradation(core_wc_model, species_types)
        wc_lang.io.Writer().run(core_wc_model, model_filename)

    def rna_weight(self, rna):
        """ Calculate molecular weight of a RNA molecule

        Args:
            rna (:obj:`string`): RNA sequence

        Returns:
            :obj:`float`: molecular weight
        """
        for j in range(len(rna)):
            rna_weight = Enrich_Polymers.RNAMW
            base = rna[j]
            rna_weight += Enrich_Polymers.BASE_WEIGHTS[base]

        return rna_weight

    def rna_charge(self, rna):
        """ Calculate charge of a RNA molecule

        Args:
            rnas (:obj:`string`): RNA sequence

        Returns:
            :obj:`int`: charge of RNA molecule
        """
        # Charge of RNA molecule is the negative of its (nucleotide length + 1)
        return -(len(rna) + 1)

    def transcription(self, core_model, species_types):
        """
        Creates transcription reactions for the randomly generated RNA sequences

        Args:
            core_model (:obj:`Model`): wc_lang random model
            species_types (:obj:'list'): list of SpeciesType objects in core_model
        """
        nameList = []
        idList = []
        for line in open(Enrich_Polymers.CORE_POLYMER_RXN_INFO):
            # TODO: simplify
            line   = line.rstrip('\n')
            lineList = line.split(',')
            ID = lineList[0]
            name = lineList[1]
            submodel = lineList[2]
            if ID == 'Id':
                continue
            if submodel == 'Transcription':
                idList.append(ID.replace('-', '_'))
                nameList.append(name)

        # write the transcription reactions into the model in memory
        # TODO: change transSubModel to trans_submodel
        transSubModel = core_model.get_component('submodel', 'Transcription')
        index = 0

        speciesList = core_model.get_species()
        atp = self.find_species(core_model,speciesList, "ATP")
        gtp = self.find_species(core_model,speciesList, "GTP")
        utp = self.find_species(core_model,speciesList, "UTP")
        ctp = self.find_species(core_model,speciesList, "CTP")
        h2o = self.find_species(core_model,speciesList, "H2O")
        ppi = self.find_species(core_model,speciesList, "PPI")
        h = self.find_species(core_model,speciesList, "H")

        for element in species_types:
            # create RNA transcription reaction with its reactants and products
            if element.type == wc_lang.core.SpeciesTypeType.rna:
                reaction = nameList[index]
                current_id = idList[index]
                transcript = wc_lang.core.Reaction(id=current_id, name = reaction, submodel = transSubModel)
                transcript.participants.create(species = atp, coefficient=-1 * element.structure.count('A'))
                transcript.participants.create(species = gtp, coefficient=-1 * element.structure.count('G'))
                transcript.participants.create(species = utp, coefficient=-1 * element.structure.count('U'))
                transcript.participants.create(species = ctp, coefficient=-1 * element.structure.count('C'))
                transcript.participants.create(species = h2o, coefficient=-1)
                transcript.participants.create(species = h, coefficient=1)
                if self.find_species(core_model,speciesList, element.id) != None:
                    transcript.participants.create(species = self.find_species(core_model,speciesList, element.id), coefficient = 1)
                else:
                    transcript.participants.create(species = wc_lang.core.Species(species_type = element,
                                                                                 compartment = core_model.compartments.get(id='c')),
                                                                                  coefficient=1)
                transcript.participants.create(species = ppi,coefficient=element.structure.count('A')+element.structure.count('G')+element.structure.count('U')+element.structure.count('C') )
                index += 1

    def translation(self, core_model, species_types):
        """ Create RNA degradation reactions for the randomly generated RNA sequences

        Args:
            core_model (:obj:`Model`): wc_lang random model
            species_types (:obj:'list'): list of SpeciesType objects in core_model
        """
        nameList = []
        idList = []
        for line in open(Enrich_Polymers.CORE_POLYMER_RXN_INFO):
            line   = line.rstrip('\n')
            lineList = line.split(',')
            ID = lineList[0]
            name = lineList[1]
            submodel = lineList[2]
            if ID == 'Id':
                continue
            if submodel == 'Translation':
                idList.append(ID.replace('-', '_'))
                nameList.append(name)

        # write the translation reactions into the model in memory
        index = 0
        submodels = core_model.get_submodels()
        for submodel in submodels: # finding the translation submodel among all the submodels
            if submodel.name == "Translation":
                transSubModel = submodel


        # getting Species objects for the common reactants and products beforehand
        speciesList = core_model.get_species()
        gtp = self.find_species(core_model,speciesList, "GTP")
        h2o = self.find_species(core_model,speciesList, "H2O")
        gdp = self.find_species(core_model,speciesList, "GDP")
        pi = self.find_species(core_model,speciesList, "PI")
        h = self.find_species(core_model,speciesList, "H")

        for aa in Enrich_Polymers.ABBREVIATIONS:  # creating Species objects for each amino acid beforehand
            if self.find_species(core_model,speciesList, Enrich_Polymers.ABBREVIATIONS[aa]) == None:
                species = wc_lang.core.Species(species_type = core_model.species_types.get(id=Enrich_Polymers.ABBREVIATIONS[aa]),
                                                                      compartment = core_model.compartments.get(id='c'))


        for element in species_types:
            if element.type == wc_lang.core.SpeciesTypeType.protein: # if SpeciesType object is protein, then create translation reaction with correct numbers of reactants and products
                # print(element)
                reaction = nameList[index]
                current_id = idList[index]
                transcript = wc_lang.core.Reaction(id=current_id, name = reaction, submodel = transSubModel)

                length = len(element.structure)

                for aa in Enrich_Polymers.ABBREVIATIONS: # goes through all amino acids, counts number of that amino acid in protein sequence, adds Species object of that amino acid to the reaction
                    count = element.structure.count(aa)
                    if count != 0:
                        transcript.participants.create(species = self.find_species(core_model,speciesList, Enrich_Polymers.ABBREVIATIONS[aa]), coefficient = -1 * count)


                transcript.participants.create(species = gtp,  coefficient=-1 * (2 * length + 3)  )
                transcript.participants.create(species = h2o, coefficient=-1 * (length + 4))
                transcript.participants.create(species = gdp, coefficient= 2 * length + 3)
                transcript.participants.create(species = pi, coefficient = 2*length + 3)
                transcript.participants.create(species = h, coefficient = 2 * length + 3)
                if self.find_species(core_model,speciesList, element.id) != None:
                    transcript.participants.create(species = self.find_species(core_model,speciesList, element.id), coefficient = 1)
                else:
                    transcript.participants.create(species = wc_lang.core.Species(species_type = element,
                                                                                 compartment = core_model.compartments.get(id='c')),
                                                                                  coefficient=1)

                index += 1

    def rna_degradation(self, core_model, species_types):
        """ Creates RNA degradation reactions for the randomly generated RNA sequences

        Args:
            core_model (:obj:`Model`): wc_lang random model
            species_types (:obj:'list'): list of SpeciesType objects in core_model
        """
        nameList = []
        idList = []
        for line in open(Enrich_Polymers.CORE_POLYMER_RXN_INFO):
            line   = line.rstrip('\n')
            lineList = line.split(',')
            ID = lineList[0]
            name = lineList[1]
            submodel = lineList[2]
            if ID == 'Id':
                continue
            if submodel == 'RnaDegradation':
                idList.append(ID.replace('-', '_'))
                nameList.append(name)

        # write the rna degradation reactions into the model in memory
        index = 0
        submodels = core_model.get_submodels()
        for submodel in submodels: # finding the rna degradation submodel among all the submodels
            if submodel.name == "RNA degradation":
                transSubModel = submodel

        # findi the Species objects for the common reactants and products beforehand
        speciesList = core_model.get_species()

        gmp = self.find_species(core_model,speciesList, "GMP")
        h2o = self.find_species(core_model,speciesList, "H2O")
        amp = self.find_species(core_model,speciesList, "AMP")
        cmp = self.find_species(core_model,speciesList, "CMP")
        ump = self.find_species(core_model,speciesList, "UMP")
        h = self.find_species(core_model,speciesList, "H")

        for element in species_types:
            if element.type == wc_lang.core.SpeciesTypeType.rna: # if SpeciesType object is RNA, then create RNA degradation reaction with correct numbers of reactants and products
                reaction = nameList[index]
                current_id = idList[index]
                transcript = wc_lang.core.Reaction(id=current_id, name = reaction, submodel = transSubModel)

                transcript.participants.create(species = gmp,  coefficient= element.structure.count('G') )
                transcript.participants.create(species = amp, coefficient= element.structure.count('A')  )
                transcript.participants.create(species = cmp, coefficient= element.structure.count('C')  )
                transcript.participants.create(species = ump, coefficient = element.structure.count('U') )
                transcript.participants.create(species = h2o , coefficient = -1 * ( element.structure.count('A')+element.structure.count('G')+element.structure.count('U')+element.structure.count('C') - 1))
                transcript.participants.create(species = h , coefficient = element.structure.count('A')+element.structure.count('G')+element.structure.count('U')+element.structure.count('C') - 1)
                if self.find_species(core_model,speciesList, element.id) != None:
                    transcript.participants.create(species = self.find_species(core_model,speciesList, element.id), coefficient = 1)
                else:
                    transcript.participants.create(species = wc_lang.core.Species(species_type = element,
                                                                                 compartment = core_model.compartments.get(id='c')),
                                                                                  coefficient=1)
                index += 1

    def find_species(self, core_model, speciesList, ID):
        '''Finds species with given SpeciesType ID

        Args:
            core_model (:obj:`Model`): wc_lang random model
            speciesList (:obj:'list'): list of Species objects in core_model

        Returns:
            :obj:`Species`: Species object corresponding to the given SpeciesType ID
        '''
        for species in speciesList:
            if species.species_type.id == ID and species.compartment == core_model.compartments.get(id='c'):
                return species

        return None
