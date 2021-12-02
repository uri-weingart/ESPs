Specific Peptides
=============

# General

This is a utility package to allow for deducing protein classification from its given amino acid sequence. 
The methodology is described in our [article](https://github.com/uri-weingart/ESPs/blob/main/Specific_Peptides_Perspective_of_Proteins.pdf/).
 It can be applied to **Enzyme Classification** (using ESPs), and classification of **GPCR** proteins (using GSPs) and **Z finger proteins** (using ZSPs).

## Usage

To use our algorithm, please clone it from Github as follows:
 
git clone https://github.com/uri-weingart/ESPs.git
 
Please install the required Python dependencies, outlined in the file requirements.txt
 
We provide a fasta test file  in the data directory named  test.fasta  containing 20 enzymes downloaded from Swissprot.

## ESP searches

Go to the **src** directory

### python SPs.py -i <your_fasta_file.fasta>
  
## GSP   searches
 
### python SPs.py -i <your_fasta_file.fasta> -d  dGPCRs.json

## ZSP   searches

### python SPs.py -i <your_fasta_file.fasta> -d  dZFs.json
 