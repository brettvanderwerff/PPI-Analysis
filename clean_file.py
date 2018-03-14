'''Script cleans the dataframe returned by the swiss_gene_name_conv script to remove non-human interactions,
chemicals, RNAs and redundant entries.
'''

import id_converter
import id_parser
import numpy as np
import pandas as pd
import swiss_gene_name_conv

def remove_redundant_entries(df):
    '''Function will leverage the IDs in the 'Interaction identifier(s)' column to detect entries that
    are redundant across databases. Function will replace rows that contain
    redundant information with NaN.
    '''
    interaction_id_list = []
    for item in df['Interaction identifier(s)']:
        interaction_id_list.append(item.upper().split('|'))
    new_list = []
    for counter, list_of_interest in enumerate(interaction_id_list):
        for nested_list in interaction_id_list[(counter + 1):]:
            if set(list_of_interest).issubset(nested_list):
                if len(list_of_interest) < len(nested_list) or len(list_of_interest) == len(nested_list):
                    interaction_id_list[counter] = np.nan
    sep = '|'
    for list in interaction_id_list:
        try:
            new_list.append(sep.join(list))
        except TypeError:
            new_list.append(list)
    return pd.concat([pd.DataFrame(new_list, columns = ['Interaction identifier(s)']),
                      df.drop(columns=['Interaction identifier(s)'])], axis=1)

def replace_pubmed(df):
    '''Replaces values in the 'Publication Identifier(s)' column of the dataframe that contain IDs other than
    pubmed IDs.
    '''
    row_list = []
    pubmed_list = []
    for item in df['Publication Identifier(s)']:
        row_list.append(item.upper().split('|'))
    for nested_list in row_list:
        new_nested_list = []
        for item in nested_list:
            if 'PUBMED' in item:
                new_nested_list.append(item)
        sep = '|'
        pubmed_list.append(sep.join(new_nested_list))
    return pd.concat([pd.DataFrame(pubmed_list, columns=['Publication Identifier(s)']),
                                   df.drop(columns=['Publication Identifier(s)'])], axis=1)

def publication_compare(df):
    '''Counts the publications in the 'Publication Identifier(s)' and 'Publication 1st author(s)'] columns.
    if there is an inconsistency in number, function will replace the 'Publication Identifier(s)' value with NaN.
    '''
    new_list = []
    pubmed_list = []
    for item in df['Publication Identifier(s)']:
        pubmed_list.append(item.upper().split('|'))
    auth_row_list = []
    for item in df['Publication 1st author(s)']:
        auth_row_list.append(item.upper().split('|'))
    for c, (p,a) in enumerate(zip(pubmed_list, auth_row_list)):
        if len(p) != len(a):
            pubmed_list[c] = np.nan
    sep = '|'
    for item in pubmed_list:
        try:
            new_list.append(sep.join(item))
        except TypeError:
            new_list.append(item)
    return pd.concat([pd.DataFrame(new_list, columns=['Publication Identifier(s)']),
                      df.drop(columns=['Publication Identifier(s)'])], axis=1)

def remove_empty(df):
    '''Will remove entire columns that contain an empty gene name or any other empty value (NaN)
     in it's rows. Empty gene names are produced from
    non-human interactions, chemicals, RNAs, or anything else that does not map to a uniprot ID.
    '''
    return df.dropna(axis=0).reset_index(drop=True)

def run(df):
    '''Calls all functions in script in order.
    '''
    redun_removed = remove_redundant_entries(df=df)
    pubmed_replaced = replace_pubmed(df=redun_removed)
    compared_pub = publication_compare(df=pubmed_replaced)
    return remove_empty(df=compared_pub)

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(df=id_parsed_df)
    gene_name_conv_df = swiss_gene_name_conv.run(df=id_converted_df)
    print(run(df=gene_name_conv_df))



