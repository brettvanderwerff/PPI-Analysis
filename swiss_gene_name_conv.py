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
    reader = pd.read_csv('HUMAN_9606_idmapping.dat', delimiter="\t", names=['UniProtKB_ID', 'Identifier', 'Value'])
    gene_name_df = reader.loc[reader['Identifier'] == 'Gene_Name']
    gene_name_df_dd = gene_name_df.drop_duplicates(subset='UniProtKB_ID', keep=False)
    # drop_duplicates ensures that only UniProtKB IDs that correspond to one common gene name are used
    id_column = df[[label + ' ID']]
    id_df = id_column.rename(columns={(label + ' ID'): 'UniProtKB_ID'})
    merged = pd.merge(id_df, gene_name_df_dd[['UniProtKB_ID', 'Value']], on='UniProtKB_ID', how='left')
    df[label + ' gene name'] = np.where(merged['Value'].notnull(), merged['Value'], np.nan)
    return df

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
    '''Calls the find_gene_name function once to create a dataframe column of common gene names of protein A and again
     to create a column for the common gene names of protein B. An 'Interactor name' column is then generated for the
     dataframe with the interactor_column function. NaN values are dropped from the returned function, which will
      mostly be interactors that are not proteins or not human proteins.
    '''
    gene_name_a_df = find_gene_name(df=df, label='Parsed A')
    gene_name_b_df = find_gene_name(df=gene_name_a_df, label='Parsed B')
    return interactor_column(df=gene_name_b_df, query_gene_name=query_gene_name).dropna(axis=0).reset_index(drop=True)

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(df=id_parsed_df)
    print(run(df=id_converted_df, query_gene_name='MST1R'))


