'''Script cleans the dataframe returned by the swiss_gene_name_conv script to remove unspecified methods, unassigned
pubmed entries, redundant entries etc.
'''
import id_converter
import id_parser
import numpy as np
import pandas as pd
import swiss_gene_name_conv

def pool_methods_apply(list):
    '''Argument for the series method 'apply' in the pool_methods function. This function groups similar methods
    together under 'umbrella' terms. Adds NaN for interactions gotten from unspecified and unknown methods.
    '''
    method_dict = pd.read_csv('method_dict.txt', delimiter='\t', dtype=str)
    method_list = []
    for item in list:
        if 'unknown' in item:
            method_list = [np.nan]
            break
        method = method_dict.loc[method_dict['Code'] == str(item[8:15])]['ID'].iloc[0]
        #there is probably a newer method that is not in the dictionary if an error is traced to this above line.
        #this is easily fixed by printing the above method and catching the methods that throw the error and updating
        #the method_dict.txt file to include those error throwing methods
        if method == 'UNSPM':
            method_list = [np.nan]
            break
        if method not in method_list:
            method_list.append(method)
    return '|'.join(method_list) if len(method_list) > 1 else method_list[0]

def pool_methods(df):
    '''Applies the pool_methods_apply function to the 'Interaction detection method(s)' column of the dataframe.
    Deletes interactions gotten from unspecified and unknown methods.
    '''
    df['Pooled Methods'] = df['Interaction detection method(s)'].str.split('|').apply(pool_methods_apply)
    return df.dropna(axis=0).reset_index(drop=True)

def filter_pubmed_ID_apply(list):
    '''Function filters the pubmed IDs present in the 'Publication Identifier(s)' column to save only pubmed ID's
    that are unambiguously assigned. Places NaN in columns that have unassigned pubmed IDs.
    This function is meant to be used via the apply method of a series.
    '''
    pubmed_only = []
    for item in list:
        if 'PUBMED:UNASSIGNED' in item.upper():
            pubmed_only = [np.nan]
            break
        elif 'PUBMED' in item.upper():
            pubmed_only.append(item.upper())
    return ('|').join(set(pubmed_only)) if len(pubmed_only) > 1 else pubmed_only[0]

def filter_pubmed_ID(df):
    '''Applies the filter_pubmed_apply function to the 'Publication Identifier(s)' column of the dataframe. Removes any
    rows that have unassigned pubmed IDs.
    '''
    df['Publication Identifier(s)'] = df['Publication Identifier(s)'].str.split('|').apply(filter_pubmed_ID_apply)
    return df.dropna(axis=0).reset_index(drop=True)

def publication_compare(df):
    '''Function compares the 'Publication Identifier(s)' and 'Publication 1st author(s)' columns of the dataframe and
    removes any rows that have an ambiguous number of publications between these columns.
    '''
    df['Publication Identifier(s)'] = np.where(df['Publication Identifier(s)'].str.split('|').apply(lambda x: len(x))
                                               != df['Publication 1st author(s)'].str.split('|').apply(lambda x: len(x))
                                               , np.nan, df['Publication Identifier(s)'])
    return df.dropna(axis=0).reset_index(drop=True)

def eliminate_duplicate_genes(df):
    '''Function will group together columns that describe the same gene(i.e. multiple databases have entries describing
    the same gene.)
    '''
    df_sort = df.groupby('Interactor name').agg({'Publication Identifier(s)': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Interaction identifier(s)': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   '#ID(s) interactor A': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'ID(s) interactor B': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Alt. ID(s) interactor A': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Alt. ID(s) interactor B': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Alias(es) interactor A': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Alias(es) interactor B': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Interaction detection method(s)': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Publication 1st author(s)': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Taxid interactor A': 'first',
                                                   'Taxid interactor B': 'first',
                                                   'Interaction type(s)': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Source database(s)': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Confidence value(s)': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')])),
                                                   'Parsed A ID': 'first',
                                                   'Parsed B ID': 'first',
                                                   'Parsed A gene name': 'first',
                                                   'Parsed B gene name': 'first',
                                                 'Pooled Methods': lambda x : '|'.join(set([z.upper() for y in x for z in y.split('|')]))
                                                 }).reset_index()
    return df_sort

def run(df):
    '''Calls all functions in script in order.
    '''
    print('Pooling methods...')
    pooled_methods = pool_methods(df=df)
    print('Methods pooled')
    print('Filtering Pubmed entries...')
    join_pubmed = filter_pubmed_ID(df=pooled_methods)
    pub_compared = publication_compare(df=join_pubmed)
    print('Pubmed entries filtered')
    print('Removing redundant entries...')
    dup_genes_eliminated = eliminate_duplicate_genes(df=pub_compared)
    print('Redundant entries removed')
    return dup_genes_eliminated

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(df=id_parsed_df)
    gene_name_conv_df = swiss_gene_name_conv.run(df=id_converted_df, query_gene_name='MST1R')
    print(run(df=gene_name_conv_df))
