'''
Script exchanges Biogrid IDs for swissprot IDs in a dataframe returned from the id_parser script.
'''
import id_parser
import pandas as pd

def replace_biogrid_id(df, label):
    '''Will replace each biogrid ID with a swissprot ID for protein A or protein B in the dataframe gotten from the
     id_parser script.
     '''
    biogrid_swiss_conversion = pd.read_csv('biogrid_swiss_conversion.csv', delimiter="\t", dtype=str)
    for item in df[label]:
        if len(biogrid_swiss_conversion.loc[biogrid_swiss_conversion['BIOGRID_ID'] == item]) == 1:
            new_id = biogrid_swiss_conversion.loc[biogrid_swiss_conversion['BIOGRID_ID'] == item,
                                                  'IDENTIFIER_VALUE'].iloc[0]
            df[label] = df[label].replace([item], new_id)
    return df

def run(df):
    '''Calls the replace biogrid_id function twice, once to replace the biogrid IDs for protein A and then again
    to replace the biogrid IDs for protein B.
    '''
    update_df_A = replace_biogrid_id(df=df, label='Parsed A ID')
    update_df_B = replace_biogrid_id(df=update_df_A, label='Parsed B ID')
    return update_df_B

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    print(run(df=id_parsed_df))











