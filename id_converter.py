'''
Script exchanges Biogrid IDs for UniProtKB IDs in a dataframe returned from the id_parser script.
'''
import id_parser
import numpy as np
import pandas as pd

def replace_biogrid_id(df,label):
    '''Will replace each Biogrid ID for a UniProtKB ID for protein A or protein B in the dataframe gotten from the
    id_parser script.
    '''
    biogrid_uniprot_conversion = pd.read_csv('ID_conversion_files/biogrid_uniprot_conversion.csv', delimiter="\t", dtype=str)
    biogrid_uniprot_conversion.drop_duplicates(subset='IDENTIFIER_VALUE', keep=False, inplace=True)
    #drop_duplicates ensures that only Biogrid IDs that correspond to one UniProtKB ID are used
    id_column = df[[label]]
    id_df = id_column.rename(columns={label: 'BIOGRID_ID'})
    merged = pd.merge(id_df, biogrid_uniprot_conversion, on='BIOGRID_ID', how='left')
    df[label] = np.where(merged['IDENTIFIER_VALUE'].notnull(), merged['IDENTIFIER_VALUE'], df[label])
    return df

def run(df):
    '''Calls the replace biogrid_id function twice, once to replace the biogrid IDs for protein A and then again
    to replace the biogrid IDs for protein B.
    '''
    print('Converting Biogrid IDs to UniProtKB IDs...')
    update_df_A = replace_biogrid_id(df=df, label='Parsed A ID')
    update_df_B = replace_biogrid_id(df=update_df_A, label='Parsed B ID')
    print('Biogrid IDs converted')
    return update_df_B

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    print(run(df=id_parsed_df))












