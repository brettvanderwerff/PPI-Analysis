import numpy as np
import pandas as pd
import re

'''Script is designed to organize several functions that use RegEx to detect Biogrid IDs and Swissprot IDs in
a psi_mitlab file from the IMEX protein-protein interaction database found at:
http://www.ebi.ac.uk/Tools/webservices/psicquic/view/
'''

def regex(expression, target):
    '''Basic RegEx function that takes a RegEx expression and a target string for parsing as arguments.
    '''
    id_regex = re.compile(expression)
    mo = id_regex.search(target)
    return mo.group()

def get_id(psi_mitlab, column, alias_column, protein_label):
    '''Gets swissprot, biogrid, entrez ids from a psi_mitlab file. The search prioritizes swissprot and biogrid IDs
    because unlike entrez IDs, swissprot and biogrid IDs usually reference one protein.
    '''
    swissprot_expression = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
    biogrid_expression = r'biogrid:(\d\d\d\d\d\d)'
    enterez_expression = r'entrez gene/locuslink:(\d)+'
    id_list = []
    for counter, item in enumerate(psi_mitlab[column], 0):
        try:
            id_list.append((regex(expression=biogrid_expression,
                           target=item))[8:])
        except AttributeError:
            try:
                id_list.append(regex(expression=swissprot_expression,
                               target=item))
            except AttributeError:
                try:
                    id_list.append(
                        regex(expression=swissprot_expression, target=psi_mitlab[alias_column][counter]))
                except AttributeError:
                    try:
                        id_list.append(
                            (regex(expression=enterez_expression, target=psi_mitlab[alias_column][counter]))[
                            22:])
                    except AttributeError:
                        id_list.append(np.nan)
    return pd.DataFrame(id_list, columns=[protein_label])

def combine_dfs(psi_mitlab, protein_A_ID, protein_B_ID):
    '''Combines dataframes with extracted labels for protein A and protein B (protein_A_ID and protein_b_ID) with the
    psi_mitlab dataframe.
    '''
    return pd.concat([psi_mitlab, protein_A_ID, protein_B_ID], axis=1)

def run(filename):
    '''Run function invokes all other functions in this script in the correct order. Accepts the filename of a
    psi_mitlab file as an argument.
    '''
    psi_mitlab = pd.read_csv(filename, delimiter='\t')
    protein_A_ID = get_id(psi_mitlab=psi_mitlab, column='#ID(s) interactor A',
                          alias_column='Alias(es) interactor A', protein_label='Parsed A ID')
    protein_B_ID = get_id(psi_mitlab=psi_mitlab, column='ID(s) interactor B',
                          alias_column='Alias(es) interactor B', protein_label='Parsed B ID')
    return combine_dfs(psi_mitlab=psi_mitlab, protein_A_ID=protein_A_ID,
                       protein_B_ID=protein_B_ID)

if __name__ == '__main__':
    print(run(filename='clusteredQuery_MST1R.txt'))











