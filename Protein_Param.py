from consensus_seq import consensus_Seq
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class protein_param:
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
    }

    def __init__(self):
        self.seq = consensus_Seq().raise_error()
        self.codon_table = dict(protein_param.DNA_codon_table)

    def translate_seq(self):
        codon_table = self.codon_table
        trans_seq = Seq(self.seq)
        protein = ""
        for i in range(0, len(trans_seq), 3):
            codon = trans_seq[i:i + 3]
            protein += codon_table[codon]
        return(protein)

    def param_analysis(self):
        analysis_seq = ProteinAnalysis(self.translate_seq())
        return(analysis_seq)

    def amino_acid_comp(self):
        comp_seq = self.param_analysis()
        comp = comp_seq.count_amino_acids()
        return(comp)

    #def amino_acid_percent(self):
        #percent_seq = self.param_analysis()
        #percent = percent_seq.get_amino_acids_percent()
        #print(percent)

    def amino_acid_mw(self):
        mw_seq = self.param_analysis()
        mw = mw_seq.molecular_weight()
        return(mw)

    def amino_acid_aromaticity(self):
        aromaticity_seq = self.param_analysis()
        aromaticity = aromaticity_seq.aromaticity()
        return(aromaticity)

    def amino_acid_instability(self):
        instability_seq = self.param_analysis()
        instability = instability_seq.instability_index()
        print(instability)

    def amino_acid_flexibility(self):
        flexiblity_seq = self.param_analysis()
        flexiblity = flexiblity_seq.flexibility()
        print(flexiblity)

    def amino_acid_iso(self):
        iso_seq = self.param_analysis()
        iso = iso_seq.isoelectric_point()
        print(iso)



protein_param().amino_acid_iso()