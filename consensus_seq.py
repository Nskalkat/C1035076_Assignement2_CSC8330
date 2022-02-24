import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

class consensus_Seq:

    def __int__(self):
        self.file_name = sys.argv[1]
        self.similarity_seq = float(sys.argv[2])




