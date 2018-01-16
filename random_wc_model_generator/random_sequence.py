# read and store codon translation table
file = 'translationTable.csv'

translation_table = {}

for line in open(file):
    if line[0] == '#': continue      
    line   = line.rstrip('\n')       

    fields = line.split(',')
    
    # key = each codon in translation table, value = corresponding amino acid
    translation_table[fields[0]] = fields[1]

# paramater length is desired length (in amino acids) of protein
# returns nucleotide sequence for each protein 
def protSeq(length):
    # start codon
    protein = 'AUG'

    for i in range(length-2):
        valid = False
        
        while not(valid):
            nuc1 = '{}'.format(np.random.choice(['A','U','C','G']))
            nuc2 = '{}'.format(np.random.choice(['A','U','C','G']))
            nuc3 = '{}'.format(np.random.choice(['A','U','C','G']))
            codon = nuc1 + nuc2 + nuc3
            # if codon exits in list of valid codons (translationTable)
            if codon in translation_table.keys():
                valid = True
        protein += codon
    
    # stop codon
    protein += 'UAG'

    return protein

def RNAseq(protein_sequence):

    rna_sequence = ''.join(protein_sequence)

    return rna_sequence

# complement of RNA sequence
def DNAseq(rna_sequence):
    dna_sequence = ''
    complement = {'A':'T','U':'A','C':'G','G':'C'}
    
    for nuc in rna_sequence:
        dna_sequence += complement[nuc]
    
    return dna_sequence

protein_sequence = []
# number of proteins
n = 3
# create list of random lengths (in terms of amino acids) for each of the n proteins
# could be potentially imported, with set length of genes
protLen = np.random.randint(50,1000,n)

for i in range(n):
    protein = protSeq(protLen[i])
    # list of nucleotide sequences for each protein
    protein_sequence.append(protein)

rna_sequence = RNAseq(protein_sequence)
dna_sequence = DNAseq(rna_sequence)
