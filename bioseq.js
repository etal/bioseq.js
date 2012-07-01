// Biomolecular sequence utilities
//
// Copyright (c) 2012 Eric Talevich. All rights reserved.
// Distributed under the new BSD license; see the LICENSE file.

// ENH:
//  - use Array functions to tighten up the code
//  - JSDoc
//  - figure out namespacing for browser-side use

// (function() {
// var bioseq = {};

// ---------------------------------------------------------------------
// Constant definitions

// Alphabets
const DnaBases = ['A', 'C', 'G', 'T'];
const RnaBases = ['A', 'C', 'G', 'U'];
const AminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'];
const AminoAcids3 = ['Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile',
                     'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser',
                     'Thr', 'Val', 'Trp', 'Tyr'];
const DnaComplements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'};

geneticcode = require('./geneticcode');

// ---------------------------------------------------------------------
// Functions on strings

// Count each type of letter in the seq
// ENH - accept an alphabet as an optional argument; initialize counts with it
module.exports.letterCounts = letterCounts = function (seq) {
    var counts = {};
    for (var i in seq) {
        chr = seq[i];
        counts[chr] = (counts[chr] || 0) + 1;
    }
    return counts;
}

// Get the frequency of each letter in the sequence (all sum to 1)
module.exports.letterFrequencies = letterFrequencies = function (seq) {
    var total = 0;
    var counts = letterCounts(seq);
    // sum the counts [like sum(counts.values())]
    for (var letter in counts) {
        total += counts[letter]
    }
    // divide each count by that sum
    var freqs = {};
    for (var letter in counts) {
        freqs[letter] = counts[letter] / total;
    }
    return freqs;
}

module.exports.complement = complement = function (seq) {
    // ENH: warn about non-DNA characters
    var out = [];
    for (var i in seq) {
        // Swap base pairs
        out.push(DnaComplements[seq[i]]);
    }
    return out.join('');
}

module.exports.reverseComplement = reverseComplement = function (seq) {
    // ENH: warn about non-DNA characters
    var out = [];
    // Traverse the string in reverse
    for (var i = seq.length - 1; i != 0; i--) {
        // Swap base pairs
        out.push(DnaComplements[seq[i]]);
    }
    return out.join('');
}

// Transcribe a DNA sequence to RNA
module.exports.transcribe = transcribe = function (seq) {
    var chr;
    var out = [];
    // Replace all T with U
    for (var i in seq) {
        chr = seq[i];
        if (chr == 'T') {
            out.push('U');
        } else if (chr == 'A' || chr == 'C' || chr == 'G') {
            out.push(chr);
        } else {
            // ENH - deal with ambiguity codes
            throw new Error("Sequence is not a DNA alphabet");
        }
    }
    return out.join('');
}

// Transcribe an RNA sequence back to DNA
module.exports.backTranscribe = backTranscribe = function (seq) {
    var chr;
    var out = [];
    // Replace all U with T
    for (var i in seq) {
        chr = seq[i];
        if (chr == 'U') {
            out.push('T');
        } else if (chr == 'A' || chr == 'C' || chr == 'G') {
            out.push(chr);
        } else {
            // ENH - deal with ambiguity codes
            throw new Error("Sequence is not an RNA alphabet");
        }
    }
    return out.join('');
}

// Translate a DNA or RNA sequence to protein
module.exports.translate = translate = function (seq, codonTableId) {
    // Default to the generic/universal codon table
    var codonTable = geneticcode.CodonTables[codonTableId || 1];
    if (codonTable == undefined) {
        throw new Error("Invalid codon table ID " + codonTableId);
    }
    var codon;
    var out = [];
    if (seq.indexOf('U') != -1) {
        // Codon tables are for DNA, not RNA, but we can compensate
        // XXX Biologically, it should be the other way...
        console.warn("Back-transcribing RNA to DNA for translation");
        seq = backTranscribe(seq);
    }
    for (var i = 0; i < Math.floor(seq.length/3); i++) {
        codon = seq.slice(i*3, (i+1)*3);
        out.push(codonTable[codon]);
    }
    return out.join('');
}

module.exports.translate3frames = translate3frames = function (seq, codonTableId) {
    // ENH - use _.map()
    return [translate(seq, codonTableId),
            translate(seq.slice(1), codonTableId),
            translate(seq.slice(2), codonTableId)];
}

module.exports.translate6frames = translate6frames = function (seq, codonTableId) {
    var forward = translate3frames(seq, codonTableId);
    var reverse = translate3frames(reverseComplement(seq), codonTableId);
    return forward.concat(reverse);
}

module.exports.ungap = ungap = function (seq) {
    return seq.split('-').join('');
}


// ---------------------------------------------------------------------
// The core object

function Sequence(data, id, description, alphabet, features, annot) {
    this.data = data;
    this.id = id;
    this.description = description;
    this.alphabet = alphabet;
    this.features = features;
    this.annot = annot;
}
// Convert the global functions into methods that apply to the data attribute
// ENH - automate this somehow
Sequence.prototype.letterCounts = function () {
    return letterCounts(this.data);
}
Sequence.prototype.letterFrequencies = function () {
    return letterFrequencies(this.data);
}
Sequence.prototype.complement = complement = function () {
    return complement(this.data);
}
Sequence.prototype.reverseComplement = function () {
    return reverseComplement(this.data);
}
Sequence.prototype.transcribe = function () {
    return transcribe(this.data);
}
Sequence.prototype.backTranscribe = function () {
    return backTranscribe(this.data);
}
Sequence.prototype.translate = function (codonTableId) {
    return translate(this.data, codonTableId);
}
Sequence.prototype.translate3frames = function (codonTableId) {
    return translate3frames(this.data, codonTableId);
}
Sequence.prototype.translate6frames = function (codonTableId) {
    return translate6frames(this.data, codonTableId);
}
Sequence.prototype.ungap = function () {
    return ungap(this.data);
}
module.exports.Sequence = Sequence;


// ---------------------------------------------------------------------
// I/O

// TODO - sort out the API for multiple formats, file I/O

module.exports.readFasta = readFasta = function (textblock) {
    var line, spaceloc, seqId, seqDescr, seqData = [];
    var sequences = [];
    var lines = textblock.split("\n");
    for (var i in lines) {
        line = lines[i];
        if (line[0] == '>') {
            // Save the sequence data that's already been read
            if (seqId) {
                sequences.push(new Sequence(seqData.join(''), seqId, seqDescr));
            }
            // Begin building a new sequence
            seqData = [];
            var spaceloc = line.indexOf(' ');
            if (spaceloc == -1) {
                // No description given
                seqId = line.slice(1).trim();
                seqDescr = "";
            } else {
                seqId = line.slice(1, spaceloc);
                seqDescr = line.slice(spaceloc + 1).trim();
            }
        } else if (line) {
            seqData.push(line.trim());
        }
    }
    // Save the final sequence
    sequences.push(new Sequence(seqData.join(''), seqId, seqDescr));
    return sequences;
}

module.exports.toFasta = writeFasta = function (seq) {
    result = ">" + seq.id;
    if (seq.description) {
        result += " " + seq.description;
    }
    // ENH: wrap lines at a given width
    result += "\n" + seq.data;
    return result;
}

//}).call();
//}).call(this);

// module.exports.bioseq = bioseq;

