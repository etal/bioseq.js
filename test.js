bs = require("./bioseq");

// Functions on raw strings
rawseq = "GATTACATG";

rnaseq = bs.transcribe(rawseq);
bs.backTranscribe(rnaseq);

bs.complement(rawseq);
bs.reverseComplement(rawseq);
bs.letterCounts(rawseq);
bs.letterFrequencies(rawseq);
bs.translate(rawseq);
bs.translate3frames(rawseq);
bs.translate6frames(rawseq);
bs.ungap(rawseq);

// Methods on the Sequence object
myseq = bs.parseFasta(">myseq hey hey hey\n" + rawseq)[0];

myrnaseq = myseq.transcribe();
myrnaseq.backTranscribe();

// myseq.complement();
myseq.letterCounts();
myseq.letterFrequencies();
myseq.reverseComplement();
myseq.translate();
myseq.translate3frames();
myseq.translate6frames();
myseq.ungap();

console.log(bs.toFasta(myseq));
