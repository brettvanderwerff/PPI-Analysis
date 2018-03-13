'''Script cleans the dataframe returned by the swiss_gene_name_conv script to remove non-human interactions,
chemicals, RNAs. This script will eventually also remove redundant entries in the dataframe, but is currently incomplete.
'''

import id_converter
import id_parser
import swiss_gene_name_conv

def remove_empty_genes(gene_name_conv_df):
    '''Will remove entire columns that contain an empty gene name in it's row. Empty gene names are produced from
    non-human interactions, chemicals, RNAs, or anything else that does not map to a uniprot ID.
    '''
    return gene_name_conv_df.dropna(axis=0)

def remove_redundant_entries(drop_na_df):
    '''Function will eventually leverage the IDs in the 'Interaction identifier(s)' column to entries that
    are redundant across databases. Currently this function is unfinished.
    '''
    interaction_id_list = []
    for item in drop_na_df['Interaction identifier(s)']:
        interaction_id_list.append(item.upper().split('|'))
    print(interaction_id_list)

def run(gene_name_conv_df):
    return remove_empty_genes(gene_name_conv_df)


if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(id_parsed_df=id_parsed_df)
    gene_name_conv_df = swiss_gene_name_conv.run(id_converted_df=id_converted_df)
    print(run(gene_name_conv_df=gene_name_conv_df))


