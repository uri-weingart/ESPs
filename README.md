SPs:   Enzymatic Specific Peptides
===================================

Utility package to predict the Enzymatic Classification (EC) number for a given amino acid sequence.

We follow our methodology developed and tested in the past, pointing out the
existence of **Specific Peptides (SPs)**  which are motifs of length ≥ 7 amino acids,
occurring on enzymes only.
 
You can refer to our [article](https://github.com/uri-weingart/ESPs/blob/main/Specific_Peptides_Perspective_of_Proteins.pdf/).

To use our algorithm, please clone it from Github as follows:
 
git clone https://github.com/uri-weingart/ESPs.git
 
Please install the required Python dependencies, outlined in the file requirements.txt
 
We provide a fasta test file  in the data directory named  test.fasta  containing 100 enzymes downloaded from Swissprot.

# To run:

### python ESPs.py -i <your_fasta_file.fasta>
  
To use the example file provided:

### python ESPs.py -i data/test.fasta
