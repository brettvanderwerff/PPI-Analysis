# Protein-Protein-Interaction-Network-Analysis

This program is intended to analyze protein protein interaction networks via analysis of data obtained from the PSICQUIC web app (http://www.ebi.ac.uk/Tools/webservices/psicquic/view/).
This is a Python implementation of the technique described in the 2017 paper "Weighted Protein Interaction Network Analysis of Frontotemporal Dementia" by Ferrari et al.
I have no affiliation with any of the authors. I built this program as a mechanism to learn Pandas and Python, as such this program and its outputs come with absolutely no warranty or guarantee. 

~~~~~~~~~~~~~~
Requirements: 
~~~~~~~~~~~~~~

numpy==1.13.3

pandas==0.22.0

**Download and place these files in the 'ID_conversion_files' folder in the main directory:**

- 'BIOGRID-IDENTIFIERS-3.4.158.tab.txt' gotten from: https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.4.158/

- 'HUMAN_9606_idmapping.dat' gotten from: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

**Download and place this file in the  'CRAPome files' folder in the main directory:**

- 'CRAPome database (H. sapiens) V 1.1 ( matrix format ).txt' gotten from: 'http://crapome.org/?q=Download'

A psit-mitlab file of your favorite protein gotten from searching for your favorite protein on 
http://www.ebi.ac.uk/Tools/webservices/psicquic/view/ and selecting Biogrid, Innate DB, MINT, and IntAct databases. Place this file in the main directory. This program was tailored to work with results from these databases specifically, but likely will work with any database that uses UniProtKB IDs/Biodgrid IDs to annotate protein interactions. The repo contains an example of this file labeled as 'clusteredQuery__MST1R.txt'

**Note: Names of files placed in the 'ID_conversion_files' and 'CRAPome files' folders must match those given above. Files downloaded from the above resources will need to be unzipped if they are downloaded in a zipped format.**

**Note: This program is only compatible with Python 3**

~~~~~~~~~~~~~~~~~~~~~~~
How to run this program:
~~~~~~~~~~~~~~~~~~~~~~~
 
 1. Clone repo
 2. Install all requirements. This includes both library dependencies and retrieving all the non-library files from the links indicated in the requirements above and placing them in the correct folders. 
 3. Execute the run function from app.py, using your psi-mitlab file name as an argument for filename and your query protein name as an argument for query_gene_name. (see app.py if __name__ == '__main__' conditional for an example.
4. Allow the program to run, it will take several minutes depending on the size of your psi-mitlab file. The first time the program is run it will install a gene ID conversion table to the directory. Generation of this table is a memory intensive process and will take several minutes to run, but it is a one time process that only is performed the first time the program is run.
5. Enjoy your results, which are returned as a dataframe of stringently filtered protein protein interactions from the psi-mitlab file. The docstrings in each script are detailed and reading through them will provide a better idea of what each function/script is doing and how to fix some errors that may pop up. I highly recommend reading the Ferrari et al. paper to fully understand what is happening. 

~~~~~~~~~~~~~~~~~~~~~
What each script does:
~~~~~~~~~~~~~~~~~~~~~

biogrid_uniprot_conv - Installs a conversion table to convert Biogrid IDs to UniprotKB IDs

id_parser - Uses RegEx to extract Biogrid IDs or UniprotKB IDs from each row of a psi-mitlab file

id_converter - Converts all Biogrid IDs to a corresponding UniprotKB ID

uniprot_gene_name_conv - Converts UniprotKB IDs to a common gene name

clean_file - Drops duplicates, non-human proteins, non-proteins, etc. from the psi-mitlab file

calc_weighted_score - Calculates the weighted protein protein interaction score

app - runs the program by calling all above scripts in order
