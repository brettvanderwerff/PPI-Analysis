# Protein-Protein-Interaction-Network-Analysis

This program is intended to analyze protein interaction networks via analysis of data obtained from the PSICQUIC web app (http://www.ebi.ac.uk/Tools/webservices/psicquic/view/).
This is a Python implementation of the technique described in the 2017 paper "Weighted Protein Interaction Network Analysis of Frontotemporal Dementia" by Ferrari et al.
I have no affiliation with any of the authors. I built this program as mechanism to learn Pandas, as such this program comes with absolutely no warranty or guarantees. 

~~~~~~~~~~~~~~
Requirements: 
~~~~~~~~~~~~~~

numpy==1.13.3

pandas==0.22.0

Biogrid ID mapping file 'BIOGRID-IDENTIFIERS-3.4.158.tab.txt' gotten from: https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.4.158/
Place this file in a folder called 'ID_conversion_files' in the main directory

'HUMAN_9606_idmapping.dat' downloaded from:
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
Place this file in a folder called 'ID_conversion_files' in the main directory

'CRAPome database (H. sapiens) V 1.1 ( matrix format ).txt' downloaded from:
'http://crapome.org/?q=Download'. Place this file in a folder called 'CRAPome files' in the main directory.

A psit-mitlab file of your favorite protein downloaded after searching for your favorite protein on 
http://www.ebi.ac.uk/Tools/webservices/psicquic/view/ and selecting Biogrid, Innate DB, MINT, and IntAct databases. Place this file in the main directory. The program was tailored to work with results from these databases, but likely will work with any database that uses UniProtKB IDs to annotate protein interactions. The repo contains an example of this file labeled as 'clusteredQuery__MST1R.txt'

~~~~~~~~~~~~~~~~~~~~~~~
How to run this program:
~~~~~~~~~~~~~~~~~~~~~~~
 
 1. Clone repo
 2. Install all requirements. This includes both library dependencies and retrieving all the non-library files from the links indicated in the requirements above and placing them in the correct directory. 
 3. Execute the run function from app.py, using your psi-mitlab file name as an argument for filename and you query protein name as an argument for query_gene_name. (see app.py if __name__ == '__main__' conditional for an example.
4. Allow the program to run, it will take several minutes depending on the size of your psi-mitlab file. The first time the program is run it will install a gene ID conversion table to the directory. Generation of this table is an intensive process and will take several minutes to run, but it is a one time process that only is performed upon the first time running the program.
5. Enjoy your results, which are returned as a dataframe of stringently filtered protein protein interactions from the psi-mitlab file. The docstrings are detailed and reading through them will provide a better idea of what each function/script is doing and how to fix some errors that may pop up. I highly recommend reading the Ferrari et al. paper to fully understand what is happening. 
