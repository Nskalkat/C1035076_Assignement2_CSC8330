from consensus_seq import consensus_Seq
from Bio.Seq import Seq

class protein_param:
    DNA_codon_table = {

    # U             C             A             G
    #U
    'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGU': 'Cys',
    'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys',
    'TTA': 'Leu', 'TCA': 'Ser', 'TAA': 'Stop', 'TGA': 'Stop',
    'TTG': 'Leu', 'TCG': 'Ser', 'TAG': 'Stop', 'TGG': 'Trp',
    #C
    'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg',
    'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
    'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
    'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
    #A
    'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser',
    'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
    'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
    'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
    #G
    'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly',
    'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
    'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
    'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
    }

    def __init__(self):
        self.seq = consensus_Seq().raise_error()

    def translate_seq(self):
        trans_seq = Seq(self.seq)
        final_seq = trans_seq.translate()
        print(final_seq)

protein_param().translate_seq()