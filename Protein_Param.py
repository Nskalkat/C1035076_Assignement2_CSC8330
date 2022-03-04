from consensus_seq import consensus_Seq #imports the consensus sequence code from previous file
from Bio.SeqUtils.ProtParam import ProteinAnalysis #code used for analysis of proteins such codon composition, PI, etc.

class protein_param:
    ''' The purpose of this class is to take the consensus DNA sequence from consensus_seq.py translate it and
        run protein analysis code on it inorder to get more information on the created consensus sequence. This class
        will also take the all the methods ran through here and then create an output file for the user to access with
        all relevant information such the consensus sequence, the methods ran in the consensus sequence python file, and
        the methods ran an this file.

        Example:
            Input:  Consensus Sequence: ATGCGTAAAGGGGAAGAACTGTTTACCGGCGTTGTGCCGATTCTGGTCGAACTGGATGGGGATGTGAAGGGGCATAAATTCTCCGTGCGTGGCGAAGGTGAAGGCGATGCAACAAACGGAAAACTGACGCTGAAATTTATCTGTACGACGGGGAAACTGCCGGTTCCGTGGCCGACGCTGGTGACCACGCTGACGTACGGCGTGCAGTGCTTTGCCCGCTATCCGGATCATATGAAGCAGCATGATTTCTTCAAAGGCGCCATGCCGGAAGGCTATGTGCAGGAGAGAACCATCTCTTTTAAAGATGATGGCACTTACAAAACCCGGGCTGAAGTGAAATTTGAAGGTGATACCCTGGTTAATCGTATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAACATTCTCGGCCATAAATTAGAATACAACTTTAATGGGCATAACGTGTATATCACGGCTGATAAACAAAAAAACGGGATTAAAGCGAACTTTAAAATTCGTCAGAATGTTGAAGAGGGCGGCGTGCAGCTTGCTGATCATTATCAGCAAAATACGCCGATTGGCGATGGTCCGGTGCTGCTGCCGGATAATCATTATCTGGGTACGCAATCTGTGCTGTCAAAAGACCCGAATGAAAAACGCGATCACATGGTGCTGCTGGAATTTGTTACGGCGGCAGGCATCACTCACGGCATGGATGAA

            Output:
            Protein Sequence: MRKGEELFTGVVPILVELDGDVKGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFARYPDHMKQHDFFKNAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNKHNVYITADKQKNGIKANFKIRQNVEEGNVQLADHYQQNTPIGDGPVLLPDNHYLNTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK


            AA composition: {'A': 9, 'C': 2, 'D': 19, 'E': 15, 'F': 16, 'G': 22, 'H': 10, 'I': 11, 'K': 20, 'L': 20, 'M': 5, 'N': 13, 'P': 10, 'Q': 8, 'R': 8, 'S': 4, 'T': 18, 'V': 18, 'W': 1, 'Y': 9}


            Molecular Weight: 27032.297900000016


            Aromaticty: 0.09243697478991597


            Instability index: 20.35840336134453


            Isoelectric value: 5.993978691101074
            '''
    DNA_codon_table = {

    # U             C             A             G
    #U
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': 'Stop', 'TGA': 'Stop',
    'TTG': 'L', 'TCG': 'S', 'TAG': 'Stop', 'TGG': 'W',
    #C
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    #A
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    #G
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
    } #Table used to translate DNA into Proteins by codons

    def __init__(self):
        'initalizes the dictionary above as well as the consensus sequence from the consensus_seq python file'
        self.seq = consensus_Seq().raise_error()
        self.codon_table = dict(protein_param.DNA_codon_table)

    def translate_seq(self):
        'The purpose of this method is to translate the DNA sequence into proteins'
        codon_table = self.codon_table
        trans_seq = self.seq
        protein_seq = ""
        for i in range(0, len(trans_seq), 3): #code used to loop to find codons in DNA seqeunce
            codon = trans_seq[i:i + 3]
            protein_seq += codon_table[codon] #finds codon in dictionary and returns key
        return(protein_seq)
        #code above adapted from https://www.geeksforgeeks.org/dna-protein-python-3/

    def param_analysis(self):
        'The purpose of this method is to run the protein sequence in the Protein param analysis to get infomration on it'
        analysis_seq = ProteinAnalysis(self.translate_seq()) #code that will allow the use protein param methods
        return(analysis_seq)

    def amino_acid_comp(self):
        'The Purpose of this method is to get amino acid composition of the protein sequence'
        comp_seq = self.param_analysis()
        comp = comp_seq.count_amino_acids()
        return(comp)

    def amino_acid_mw(self):
        'The Purpose of this method is to get molecular weight of the protein sequence'
        mw_seq = self.param_analysis()
        mw = mw_seq.molecular_weight()
        return(mw)

    def amino_acid_aromaticity(self):
        'The Purpose of this method is to get aromaticity of the protein sequence'
        aromaticity_seq = self.param_analysis()
        aromaticity = aromaticity_seq.aromaticity()
        return(aromaticity)

    def amino_acid_instability(self):
        'The Purpose of this method is to get the instability index of the protein sequence'
        instability_seq = self.param_analysis()
        instability = instability_seq.instability_index()
        return(instability)

    def amino_acid_iso(self):
        'The Purpose of this method is to isoelectric point of the protein sequence'
        iso_seq = self.param_analysis()
        iso = iso_seq.isoelectric_point()
        return(iso)

    def write_to_file(self):
        'The purpose of this method is to create an output file that the user can easily access with all relevant info'
        with open('output.txt', 'w') as f:
            f.write('Consensus Sequence: ')
            f.writelines(''.join(consensus_Seq().raise_error())) #writes the consenus sequence
            f.write('\n\n\n')
            f.write('GC percent: ')
            f.writelines(''.join(str(consensus_Seq().count_GCcomp()))) #writes GC comp
            f.write('\n\n\n')
            f.write('Nucleotide Composition in order of G,C,A,T: ')
            f.writelines(''.join(str(consensus_Seq().nucleotide_comp()))) #writes nucleotide comp
            f.write('\n\n\n')
            f.write('Protein Sequence: ')
            f.writelines(''.join(protein_param().translate_seq()))# writes protein sequence
            f.write('\n\n\n')
            f.write('AA composition: ')
            f.writelines(''.join(str(protein_param().amino_acid_comp()))) #writes aa composition
            f.write('\n\n\n')
            f.write('Molecular Weight: ')
            f.writelines(''.join(str(protein_param().amino_acid_mw()))) #writes molecular weight
            f.write('\n\n\n')
            f.write('Aromaticty: ')
            f.writelines(''.join(str(protein_param().amino_acid_aromaticity()))) #writes aromaticity
            f.write('\n\n\n')
            f.write('Instability index: ')
            f.writelines(''.join(str(protein_param().amino_acid_instability()))) #writes instability index
            f.write('\n\n\n')
            f.write('Isoelectric value: ')
            f.writelines(''.join(str(protein_param().amino_acid_iso()))) #writes Isoelectric value




