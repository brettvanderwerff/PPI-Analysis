'''
Script locates common gene names for UniProtKB IDs in a dataframe returned from the id_converter script. The dictionary
referenced is 'HUMAN_9606_idmapping.dat' downloaded from:
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
'''
import id_converter
import id_parser
import numpy as np
import pandas as pd

def find_gene_name(df, label):
    '''Finds the corresponding common gene name for each human UniProtKB ID in the df returned by the run function
    of id_converter. If gene name is not found then it is likely a not human or not a protein.
    '''
    gene_name_reader = pd.read_csv('HUMAN_9606_idmapping.dat', delimiter="\t", header=None)
    gene_name_only_df = gene_name_reader.loc[gene_name_reader[1] == 'Gene_Name']
    gene_name_list = []
    for item in df[label + ' ID']:
        try:
            gene_name_list.append(gene_name_only_df.loc[gene_name_only_df[0] == item][2].iloc[0])
        except:
            gene_name_list.append(np.nan)
    return pd.DataFrame(gene_name_list, columns=[label + ' gene name'])

def interactor_column(df, query_gene_name):
    '''Generates a column that specifies only the interactors of the query protein. i.e. if the query protein is
    MST1R and the interactor is MAPK, the MAPK would populate the "Interactor Name' column.
    '''
    df['Interactor name'] = np.where(df['Parsed A gene name'] != query_gene_name,
                                                    df['Parsed A gene name'], np.nan)
    df['Interactor name'] = np.where(df['Parsed B gene name'] != query_gene_name,
                                                    df['Parsed B gene name'],
                                                    df['Interactor name'])
    df['Interactor name'] = np.where(
        (df['Parsed A gene name'] == query_gene_name) & (df['Parsed B gene name'] == query_gene_name),
        df['Parsed A gene name'], df['Interactor name'])
    return df

def run(df, query_gene_name):
    '''Calls the find_gene_name function once to create a dataframe of common gene names of protein A and again to
    create a dataframe for the common gene names of protein B. These dataframes are then concatenated with the dataframe
    returned by the id_converter module run function. An 'Interactor name' column is then generated for the dataframe
     with the interactor_column function. NaN values are dropped from the returned function, which will mostly be
     interactors that are not proteins or not human proteins.
    '''
    gene_name_a_df = find_gene_name(df=df, label='Parsed A')
    gene_name_b_df = find_gene_name(df=df, label='Parsed B')
    concat_dfs = pd.concat([df, gene_name_a_df, gene_name_b_df], axis=1)
    return interactor_column(df=concat_dfs, query_gene_name=query_gene_name).dropna(axis=0).reset_index(drop=True)

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(df=id_parsed_df)
    print(run(df=id_converted_df, query_gene_name='MST1R'))


