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
        protein_seq = ""
        for i in range(0, len(trans_seq), 3):
            codon = trans_seq[i:i + 3]
            protein_seq += codon_table[codon]
        return(protein_seq)

    def param_analysis(self):
        analysis_seq = ProteinAnalysis(self.translate_seq())
        return(analysis_seq)

    #def amino_acid_comp(self):
        #comp_seq = self.param_analysis()
        #comp = comp_seq.count_amino_acids()
        #print(comp)

    def amino_acid_percent(self):
        percent_seq = self.param_analysis()
        percent = percent_seq.get_amino_acids_percent()
        return(percent)

    def print_percent(self):
        print_per = self.amino_acid_percent()
        format=""
        for i in self.amino_acids_content_percent().values():
            format = f"{print_per:%}"
        print(format)

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
        return(instability)

    def amino_acid_flexibility(self):
        flexiblity_seq = self.param_analysis()
        flexiblity = flexiblity_seq.flexibility()
        return(flexiblity)

    def amino_acid_iso(self):
        iso_seq = self.param_analysis()
        iso = iso_seq.isoelectric_point()
        return(iso)

    def write_to_file(self):
        with open('output.txt', 'w') as f:
            f.write('Consensus Sequence: ')
            f.writelines(''.join(consensus_Seq().raise_error()))
            f.write('\n')
            f.write('\n')
            f.write('\n')
            f.write('GC percent: ')
            f.writelines(''.join(str(consensus_Seq().count_GCcomp())))
            f.write('\n')
            f.write('\n')
            f.write('\n')
            f.write('Nucleotide Composition in order of G,C,A,T: ')
            f.writelines(''.join(str(consensus_Seq().nucleotide_comp())))
            f.write('\n')
            f.write('\n')
            f.write('\n')
            f.write('Protein Sequence: ')
            f.writelines(''.join(protein_param().translate_seq()))
            f.write('\n')
            f.write('\n')
            f.write('\n')
            f.write('AA composition: ')
            f.writelines(''.join(protein_param().translate_seq()))



protein_param().print_percent()