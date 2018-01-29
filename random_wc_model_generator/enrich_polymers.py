"""
Takes in a wc_lang model as input, matches RNA/protein IDs with random RNA/protein sequences, molecular weights, and charges, outputs the updated model in an Excel spreadsheet

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-01-26
:Copyright: 2018, Karr Lab
:License: MIT
"""

from random_wc_model_generator.random_polymer import RandomSeqGen
import wc_lang.io
import wc_lang.core

class Enrich_Polymers (object):

    def __init__ (self, filename): #filename that contains wc_lang model
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
            
        # add other components to the model
        # reactions
        # Transcription, translation

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

