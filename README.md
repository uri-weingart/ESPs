Specific Peptides
=================

# General

This is a utility package to allow for deducing protein classification from its given amino acid sequence. 
The methodology is described in our [article](https://github.com/uri-weingart/SPs/blob/main/Specific_Peptides_Perspective_of_Proteins.pdf/).
 It can be applied to **Enzyme Classification** (using ESPs), and classification of **GPCR** proteins (using GSPs) and **Z finger proteins** (using ZSPs).

## Usage


To use our algorithm, please clone it from Github as follows:
 
git clone https://github.com/uri-weingart/SPs.git
 
Please install the required Python dependencies, outlined in the file requirements.txt

Please unzip the following files:

1. ESPs_Automaton_Pickle_Augmented.zip

2. ESPs_Augmented.zip

3. ESPs.zip


We offer two run modes:
=======================

**1. Quick** 

**2. Extended**

Quick run mode
==============

Use the following command:

**python Predict_ECs.py** 

This will run the search algorithm against the Positive Test file, which contains 25,309 enzymes annotated in SwissProt (See the article for its detailed description)

To run against another fasta file use the following command:

**python Predict_ECs.py  -i  Your_Input_Fasta_File** 


Extended run mode
=================

 
We provide provide three samples in the **src/data** directory  containing test fasta files taken from Swissprot:

**test.fasta** to demonstrate **Enzymatic** sample searches. Results for this search can be reviewed [here](https://github.com/uri-weingart/SPs/blob/main/ESPs_Hits_Summary.pdf/).

**GPCR_Test.fasta** to demonstrate **GPCR** sample searches.  [GPCR_Results](https://github.com/uri-weingart/SPs/blob/main/GPCR_Hits_Summary.pdf/).

**ZFs_Test.fasta** to demonstrate  **Zinc Fingers** sample searches. [ZF_Results](https://github.com/uri-weingart/SPs/blob/main/ZF_Hits_Summary.pdf/).

## ESP searches

### In the src directory you can issue any of the following commands:

python SPs.py -i <your_fasta_file.fasta>
  
## GSP   searches
 
python SPs.py -i <your_fasta_file.fasta> -d  dGPCRs.json

### Example:

python SPs.py -i src/data/ GPCR_Test.fasta  -d  dGPCRs.json

## ZSP   searches

python SPs.py -i <your_fasta_file.fasta> -d  dZFs.json
 
### Example:

python SPs.py -i  src/data/ZFs_Test.fasta -d  dZFs.json