ESPs:   Enzymatic Specific Peptides
===================================

Utility package to predict the Enzymatic Classification (EC) number for a given amino acid sequence.

We follow our methodology developed and tested in the past, pointing out the
existence of Specific Peptides (SPs) which are motifs of length ≥ 7 amino acids,
occurring on enzymes only.

To use our algorithm, please clone it from Github as follows:
 
 
** git clone https://github.com/uri-weingart/ESPs.git **
 
Please install the following python dependencies:
 
- sys
- pdb
- json
- Bio
- datetime
- platform
- re
- panda
- numpy
- string
- collections
- argparse
- toolz
- aoperator
- functools
- treelib
   
We provide a fasta test file named ** test.fasta** containing 100 enzymes downloaded from Swissprot.

They can be located in the data directory 
