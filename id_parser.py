'''Script is designed to  use RegEx to detect Biogrid IDs and UniProtKB IDs in
a psi_mitlab file from the IMEX protein-protein interaction database found at:
http://www.ebi.ac.uk/Tools/webservices/psicquic/view/
This program focuses on psi-mitlab files produced from interrogating Biogrid, Innate DB, MINT, and IntAct databases
with the psicquic tool.
'''
import numpy as np
import pandas as pd
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning) #Silences future warning

def get_id(df):
    '''Searches the dataframe for Biogrid IDs and UniProtKB IDs. The Innate DB, MINT, and IntAct databases seem to
     typically have UniProtKB IDs in either the 'ID(s) interactor' column or the 'Alias(es) interactor' column.
      Biogrid IDs can often be mapped to a single UniProtKB ID using a dictionary from the Biogrid website. Biogrid IDs
      will eventually be converted to UniProtKB IDs using the id_converter script.
      '''
    df['Parsed A ID'] = df['#ID(s) interactor A'].str.extract(r'biogrid:(\d{6})')
    df['Parsed A ID'] = np.where(
        df['Parsed A ID'].isnull(),
        df['#ID(s) interactor A'].str.extract('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})')[
            0],
        df['Parsed A ID'])
    df['Parsed A ID'] = np.where(df['Parsed A ID'].isnull(),
                                 df['Alias(es) interactor A'].str.extract(
                                     '([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})')[0],
                                 df['Parsed A ID'])
    df['Parsed B ID'] = df['ID(s) interactor B'].str.extract(r'biogrid:(\d{6})')
    df['Parsed B ID'] = np.where(df['Parsed B ID'].isnull(),
        df['ID(s) interactor B'].str.extract('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})')[
            0],
        df['Parsed B ID'])
    df['Parsed B ID'] = np.where(df['Parsed B ID'].isnull(),
                                 df['Alias(es) interactor B'].str.extract(
                                     '([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})')[0],
                                 df['Parsed B ID'])
    return df

def run(filename):
    '''Runs the script.
    '''
    print('Parsing IDs...')
    df = pd.read_csv(filename, delimiter='\t')
    parsed_df = get_id(df=df)
    print('IDs parsed')
    return parsed_df

if __name__ == '__main__':
    print(run(filename='clusteredQuery_MST1R.txt'))













