'''
Script locates common gene names for swissprot IDs in a dataframe returned from the id_converter script. The dictionary
referenced is 'HUMAN_9606_idmapping.dat' downloaded from:
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
'''
import id_converter
import id_parser
import numpy as np
import pandas as pd

def find_gene_name(df, label):
    '''Finds the corresponding common gene name for each human swissprot ID in the df returned by the run function
    of id_converter.
    '''
    gene_name_reader = pd.read_csv('HUMAN_9606_idmapping.dat', delimiter="\t", header=None)
    gene_name_only_df = gene_name_reader.loc[gene_name_reader[1] == 'Gene_Name']
    gene_name_list = []
    for item in df[label]:
        try:
            gene_name_list.append(gene_name_only_df.loc[gene_name_only_df[0] == item][2].iloc[0])
        except:
            gene_name_list.append(np.nan)
    return pd.DataFrame(gene_name_list, columns=[label + ' gene name'])

def run(df):
    '''Calls the find_gene_name function once to convert protein A swissprot ID to common gene names and again to
     convert protein B swissprot ID to common gene names and returns a concatenated dataframe with the newly added
    gene names.
    '''
    gene_names_A = find_gene_name(df=df, label='Parsed A ID')
    gene_names_B = find_gene_name(df=df, label='Parsed B ID')
    return pd.concat([df, gene_names_A, gene_names_B], axis=1)

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(df=id_parsed_df)
    run(df=id_converted_df)
