'''Script calculates the weighted protein-protein interaction score.
'''
import clean_file
import id_converter
import id_parser
import numpy as np
import pandas as pd
import swiss_gene_name_conv

def publication_score(df):
    '''Counts the number of publications in the 'Publication Identifier(s)' column to get the 'publication score'
    '''
    df['Publication Score'] = df['Publication Identifier(s)'].str.split('|').apply(lambda x: len(x))
    return df

def method_score(df):
    '''Counts the number of methods in the 'Pooled Methods' column to get the 'methods score'
    '''
    df['Method Score'] = df['Pooled Methods'].str.split('|').apply(lambda x: len(x))
    return df

def crapome_score_apply(row):
    '''Argument for the series 'apply' method in the crapome_score function. Gives the 'crapome score' based on how
    often the protein shows up in APMS experiments as determined by records downloaded from
    'http://crapome.org/?q=Download' in the form of the file:
     'CRAPome database (H. sapiens) V 1.1 ( matrix format ).txt'
     Also factors in how many different methods besides APMS were used to detect the
     interaction.
    '''
    crapome_df = pd.read_csv('CRAPome files/CRAPome database (H. sapiens) V 1.1 ( matrix format ).txt',
                               delimiter='\t')
    matrix = crapome_df.loc[crapome_df['Gene'] == row['Interactor name']].as_matrix()
    method_list =[]
    try:
        method_list.append(row['Pooled Methods'].str.split('|'))
    except:
        method_list.append(row['Pooled Methods'])
    try:
        crapome_ratio = np.count_nonzero(matrix[0][3:]) / len(matrix[0][3:])
    except IndexError:
        crapome_ratio = 0
    if 'APMS' in method_list and len(method_list) == 1 and crapome_ratio > .5:
        crapome_score = -1
    elif 'APMS' in method_list and len(method_list) == 2 and crapome_ratio > .5:
        crapome_score = -.5
    elif 'APMS' in method_list and len(method_list) > 3 and crapome_ratio > .5:
        crapome_score = 0
    elif 'APMS' in method_list and len(method_list) ==1 and crapome_ratio > .3 and crapome_ratio <.5:
        crapome_score = -.5
    else:
        crapome_score = 0
    return crapome_score

def crapome_score(df):
    '''Calculates the 'crapome score' by calling the crapome_score_apply function.
    '''
    df['Crapome Score'] = df.apply(crapome_score_apply, axis=1)
    return df

def weighted_score(df):
    '''Calculates weighted protein protein interaction score by summing the publication, method,
    and crapome scores.
    '''
    df['Weighted_score'] = df['Publication Score'] + df['Method Score'] + df['Crapome Score']
    return df

def write_results(df):
    '''Removes entries from dataframe that have a weighted score of 2 or less as determined by the weighted_score
    function.
    '''
    df['Publication Identifier(s)'] = np.where((df['Weighted_score'] > 2), df['Publication Identifier(s)'], np.nan)
    return df.dropna(axis=0).reset_index(drop=True)

def run(df):
    '''Calls all functions in script in order.
    '''
    print('Calculating weighted score...')
    scored_pub = publication_score(df=df)
    scored_method = method_score(df=scored_pub)
    scored_crapome = crapome_score(df=scored_method)
    final_score = weighted_score(df=scored_crapome)
    print('Weighted score calculated, now writing results')
    return write_results(df=final_score)

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(df=id_parsed_df)
    gene_name_conv_df = swiss_gene_name_conv.run(df=id_converted_df, query_gene_name='MST1R')
    cleaned_file_df = clean_file.run(df=gene_name_conv_df)
    print(run(df=cleaned_file_df))



