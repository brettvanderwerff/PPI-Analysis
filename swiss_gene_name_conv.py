import id_converter
import id_parser
import pandas as pd


def find_gene_name(id_converted_df, label):
    '''Finds the corresponding common gene name for each human swissprot ID in the df returned by the run function
    of id_converter.
    '''
    gene_name_reader = pd.read_csv('HUMAN_9606_idmapping.dat', delimiter="\t", header=None)
    gene_name_only_df = gene_name_reader.loc[gene_name_reader[1] == 'Gene_Name']
    for item in id_converted_df[label]:
        print(gene_name_only_df.loc[gene_name_only_df[0] == item])

if __name__ == '__main__':
    id_parsed_df = id_parser.run(filename='clusteredQuery_MST1R.txt')
    id_converted_df = id_converter.run(id_parsed_df=id_parsed_df)
    find_gene_name(id_converted_df=id_converted_df, label='Parsed A ID')