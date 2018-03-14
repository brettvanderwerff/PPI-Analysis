'''Script calculates the weighted protein-protein interaction score.
'''

import clean_file
import id_converter
import id_parser
import pandas as pd
import swiss_gene_name_conv

def publication_score(df):
    '''Calculates the publication score by counting how many publications are attributed to each interaction.
    '''
    auth_row_list = []
    ps = []
    for item in df['Publication 1st author(s)']: #might want to convert list to sets
        auth_row_list.append(item.upper().split('|'))
    for item in auth_row_list:
        ps.append(len(item))
    return pd.concat([pd.DataFrame(ps, columns=['Publication Score']), df], axis=1)

def pool_methods(df):
    '''Pools together method codes that reflect highly similar methods., results from this column are used to calculate
    the method score.
    '''
    method_dict = pd.read_csv('method_dict.txt', delimiter='\t', dtype=str)
    new_list = []
    method_list = []
    pooled_list = []
    for item in df['Interaction detection method(s)']:
        method_list.append(item.upper().split('|'))
    for item in method_list:
        row_list = []
        for element in item:
            method = method_dict.loc[method_dict['Code'] == str(element[8:15])]['ID'].iloc[0]
            if method not in row_list:
                row_list.append(method)
        new_list.append(row_list)
    sep = '|'
    for item in new_list:
        try:
            pooled_list.append(sep.join(item))
        except TypeError:
            pooled_list.append(item)
    return pd.concat([pd.DataFrame(pooled_list, columns=['Pooled Methods']), df], axis=1) # need to handle unspecified methods somehow

def method_score(df):
    '''Calculates the method score by counting how many different pooled methods were used for each interaction.
    '''
    method_list = []
    count = []
    for item in df['Pooled Methods']:
        method_list.append(item.upper().split('|'))
    for item in method_list:
        count.append(len(item))
    return pd.concat([pd.DataFrame(count, columns=['Method Score']), df], axis=1)

def run(df):
    '''Calls all functions in script in order.
    '''
    scored_pub = publication_score(df)
    pooled_meth = pool_methods(scored_pub)
    return method_score(pooled_meth)

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(df=id_parsed_df)
    gene_name_conv_df = swiss_gene_name_conv.run(df=id_converted_df)
    cleaned_file_df = clean_file.run(df=gene_name_conv_df)
    print(run(df=cleaned_file_df))
