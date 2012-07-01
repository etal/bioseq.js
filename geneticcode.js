const CodonTables = {
    1: {
        'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
        'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
        'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
        'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F',
        // register_ncbi_table(name = 'Standard',
        //                     alt_name = 'SGC0', id = 1,
        //                     stop_codons = [ 'TAA', 'TAG', 'TGA', ],
        //                     start_codons = [ 'TTG', 'CTG', 'ATG', ]
    },
    // TODO - define all the genetic codes from NCBI
};
module.exports.CodonTables = CodonTables;

