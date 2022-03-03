import sys
import random
from Bio import AlignIO
from Bio.Align import AlignInfo

class consensus_Seq:
    ''' The purpsoe of this class is to first read a file containing a number
        of fasta DNA sequences, and then from there returning a consensus sequence
        as well as its GC content and Nucleotide Composition

        Example:
            INPUT:
            Sequence 1: AATTG
            Sequence 2: ATACG
            Sequence 3: TATTC

            OUTPUT:
            Consensus Sequence = AATTG
            GC content: 0.2
            Nucleotide Composition in order of G,C,A,T: 1,0,2,2

            '''

    def make_consensus_seq(self):
        'The purpose of this method is to make the consensus sequence from a file of DNA sequenes in FASTA format'
        alignment = AlignIO.read(sys.argv[1], 'fasta') #code reads the contents of the FASTA file
        consensus_align = AlignInfo.SummaryInfo(alignment)
        #code above gives information about the alignement as well as aligning sequences
        consensus_seq = consensus_align.dumb_consensus(float(sys.argv[2]))
        #code above that makes the consensus sequence from a consensus value given by the user (from 0-1)
        upper_consensus = consensus_seq.upper()
        return(upper_consensus)


    def lower_similiarity(self):
        ' The purpose of this method is to then replace all the X values with N representing ambiogous nucleotides'
        simi_seq = self.make_consensus_seq() #calls for the previous method
        str_seq = str(simi_seq) #makes consensus sequence a string just in case it was not
        replaced = str_seq
        for i in str_seq: #code to run a loop to find all Xs and replace them with Ns
            if i == "X":
                replaced=str_seq.replace("X","N")
        return(replaced)

    def randomized_seq(self):
        '''The purpose of this method is to replace all Ns with a nucleotide so the other functions of this code can
        give accurate results'''
        replace_seq = self.lower_similiarity()
        random_list = ["A","T","C","G"] #make a list value with all possible nucleotides
        rando = random.choice(random_list)
        #random.choice allows the option of a random nucleotide to be selected from the list
        random_seq = replace_seq
        for i in random_seq: #code to run loop to change N with random nucleotide
            if i == "N":
                random_seq = replace_seq.replace("N",rando) #replaces N with a random nucleotide
        return(random_seq)


    def raise_error(self):
        'The purpose of this method is raise an N value error on the chance that there is a N that occurs in the sequence'
        undefined_seq = self.randomized_seq()
        for n in undefined_seq:
            if n == "N":
                raise Exception("N detected in Sequence, please lower threshold value") #gives error that N was found in seq
            else:
                return(undefined_seq)

    def count_GCcomp(self):
        'The purpose of this method is to measure GC content of the consensus sequence'
        count_seq = self.raise_error()
        count_GC = count_seq
        for i in count_seq:
            if i == "G" or "C":
                count_G = count_seq.count("G")
                count_C = count_seq.count("C")
                count_GC = (count_G + count_C)/ len(count_seq)
        return(count_GC)


    def nucleotide_comp(self):
        'The purpose of this method is to give the nucleotide composition of the consensus sequence'
        comp_seq = self.raise_error()
        comp_G = comp_seq
        comp_C = comp_seq
        comp_A = comp_seq
        comp_T = comp_seq
        for i in comp_seq:
            if i == "G" or "C" or "A" or "T":
                comp_G = comp_seq.count("G")
                comp_C = comp_seq.count("C")
                comp_A = comp_seq.count("A")
                comp_T = comp_seq.count("T")
        return(comp_G,comp_C,comp_A,comp_T)









