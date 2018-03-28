import pandas as pd
import os.path

def check_install():
    ''' Checks if 'biogrid_swiss_conversion.csv' is installed to the directory, returns True if it is.
    '''
    print('Checking if biogrid_swiss_conversion.csv is installed...')
    return True if os.path.isfile('biogrid_swiss_conversion.csv') else False

def get_biogrid_swiss_id_conv(in_filename, out_filename):
    '''Processes Biogrid ID mapping file 'BIOGRID-IDENTIFIERS-3.4.158.tab.txt' gotten from:
        https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.4.158/
        to obtain a csv for cross-referencing  biogrid and swissprot IDs. Process is very long and takes several
        minutes.
        '''
    read_file = pd.read_csv(in_filename, delimiter='\t', skiprows=28)
    read_file_trimmed = read_file.loc[(read_file['IDENTIFIER_TYPE'] == 'SWISS-PROT') &
                                      (read_file['ORGANISM_OFFICIAL_NAME'] == 'Homo sapiens')]
    biogrid_swissprot_ids = read_file_trimmed[['BIOGRID_ID' , 'IDENTIFIER_VALUE']]
    return biogrid_swissprot_ids.to_csv(out_filename, sep='\t', index=False)

def run():
    '''Calls all functions in script in order.
    '''
    if check_install() == True:
        print('Biogrid swissprot conversion file is already installed')
    else:
        print('Installing Biogrid swissprot conversion file, this will take several minutes...')
        get_biogrid_swiss_id_conv(in_filename='BIOGRID-IDENTIFIERS-3.4.158.tab.txt',
                                  out_filename='biogrid_swiss_conversion.csv')
        print('Biogrid swissprot conversion file is now installed')

if __name__ ==  '__main__':
    run()
