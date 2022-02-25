import sys
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
        print(upper_consensus)


    def lower_similiarity(self):
        simi_seq = self.make_consensus_seq()
        str_seq = str(simi_seq)
        for i in str_seq:
            if i == "X":
                replaced=simi_seq.replace("X","N")
        return(replaced)

    #def randomized_seq(self):
       # replace_seq = self.lower_similiarity()
        #random_list = ["A","T","C","G"]
        #replace_letter = "N"
        #list_shuffle = shuffle
        #for i in replace_seq:



    #def raise_error(self):
        #undefined_seq = self.lower_similiarity()
        #for n in undefined_seq:
            #if n == "N":
               # raise Exception("N detected in Sequence, please lower threshold value")

    def count_GCcomp(self):
        count_seq = self.lower_similiarity()
        for i in count_seq:
            if i == "G" or "C":
                count_G = count_seq.count("G")
                count_C = count_seq.count("C")
                count_GC = (count_G + count_C)/ len(count_seq)
        return(count_GC)


    def nucleotide_comp(self):
        comp_seq = self.lower_similiarity()
        for i in comp_seq:
            if i == "G" or "C" or "A" or "T":
                comp_G = comp_seq.count("G")
                comp_C = comp_seq.count("C")
                comp_A = comp_seq.count("A")
                comp_T = comp_seq.count("T")
        print(comp_seq)
        print(comp_G)
        print(comp_C)
        print(comp_A)
        print(comp_T)
        print(len(comp_seq))


consensus_Seq().make_consensus_seq()





