import pandas as pd

def get_swiss_entrez_conv(in_filename, out_filename):
    '''Processes Uniprot ID mapping file 'HUMAN_9606_idmapping_selected.tab' gotten from:
    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
    to obtain a csv for cross-referencing  uniprot and entrez ids.
    '''
    read_file = pd.read_csv(('ID_conversion_files/' + in_filename), delimiter='\t', header=None, dtype=str)
    swiss_entrez = read_file[[0,2]] #only grab uniprot and entrez id columns
    swiss_entrez.columns = ['uniprot_id', 'entrez_id']
    swiss_entrez_drop_na = swiss_entrez.dropna() # drops columns missing matching uniprot or entrez id
    return swiss_entrez_drop_na.to_csv(out_filename, sep='\t', index=False)

if __name__ == '__main__':
    get_swiss_entrez_conv(in_filename='HUMAN_9606_idmapping_selected.tab', out_filename='uni_entrez_conversion.csv')







