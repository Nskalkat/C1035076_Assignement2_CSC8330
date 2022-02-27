import sys
import random
from Bio import AlignIO
from Bio.Align import AlignInfo

class consensus_Seq:

    def __int__(self):
        self.file = sys.argv[1]
        self.similarity = float(sys.argv[2])

    def make_consensus_seq(self):
        alignment = AlignIO.read(sys.argv[1], 'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus_seq = summary_align.dumb_consensus(float(sys.argv[2]))
        upper_consensus = consensus_seq.upper()
        return(upper_consensus)


    def lower_similiarity(self):
        simi_seq = self.make_consensus_seq()
        str_seq = str(simi_seq)
        replaced = str_seq
        for i in str_seq:
            if i == "X":
                replaced=str_seq.replace("X","N")
        return(replaced)

    def randomized_seq(self):
        replace_seq = self.lower_similiarity()
        random_list = ["A","T","C","G"]
        rando = random.choice(random_list)
        random_seq = replace_seq
        for i in random_seq:
            if i == "N":
                random_seq = replace_seq.replace("N",rando)
        return(random_seq)


    def raise_error(self):
        undefined_seq = self.randomized_seq()
        for n in undefined_seq:
            if n == "N":
                raise Exception("N detected in Sequence, please lower threshold value")
            else:
                return(undefined_seq)

    def count_GCcomp(self):
        count_seq = self.raise_error()
        count_GC = count_seq
        for i in count_seq:
            if i == "G" or "C":
                count_G = count_seq.count("G")
                count_C = count_seq.count("C")
                count_GC = (count_G + count_C)/ len(count_seq)
        return(count_GC)


    def nucleotide_comp(self):
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









