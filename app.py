'''Protein-Protein-Interaction-Network-Analysis main app file runs the protein protein interaction program
by calling all the needed scripts in order via the run function. The run function takes the name of a psi-mitlab
file downloaded from http://www.ebi.ac.uk/Tools/webservices/psicquic/view/ as an argument and also takes the
the common gene name of the query gene as an argument. Depending of the size of the psi mitlab-file this
program can take several minutes to run. The program will take several additional minutes to run if it is being run
for the first time.
'''
import biogrid_uniprot_conv
import calc_weighted_score
import clean_file
import id_converter
import id_parser
import uniprot_gene_name_conv

def run(filename, query_gene_name):
    biogrid_uniprot_conv.run()
    id_parsed_df = id_parser.run(filename)
    id_converted_df = id_converter.run(df=id_parsed_df)
    gene_name_conv_df = uniprot_gene_name_conv.run(df=id_converted_df, query_gene_name=query_gene_name)
    cleaned_file_df = clean_file.run(df=gene_name_conv_df)
    return calc_weighted_score.run(df=cleaned_file_df)

if __name__ == '__main__':
    print(run(filename='clusteredQuery__MST1R.txt', query_gene_name='MST1R'))

