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
        return(upper_consensus)

    def lower_similiarity(self):
        simi_seq = self.make_consensus_seq()
        str_seq = str(simi_seq)
        for i in str_seq:
            if i == "X":
                replaced=simi_seq.replace("X","N")
        return(replaced)

    #def raise_error(self):
        #undefined_seq = self.lower_similiarity()
        #for n in undefined_seq:
            #if n == "N":
                #raise Exception("N detected in Sequence, please lower threshold value")

    def count_GCcomp(self):
        count_seq = self.lower_similiarity()
        for i in count_seq:
            if i == "G" or "C":
                count_G = count_seq.count("G")
                count_C = count_seq.count("C")
                count_GC = (count_G + count_C)/ len(count_seq)
        print(count_seq)
        percent = "GC content:{:.0%}".format(count_GC)
        print(percent)


consensus_Seq().count_GCcomp()





