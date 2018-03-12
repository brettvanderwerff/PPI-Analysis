import pandas as pd

def get_biogrid_swiss_conv(in_filename, out_filename):
    '''Processes Biogrid ID mapping file 'BIOGRID-IDENTIFIERS-3.4.158.tab.txt' gotten from:
        https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.4.158/
        to obtain a csv for cross-referencing  biogrid and swissprot IDs.
        '''
    read_file = pd.read_csv(in_filename, delimiter='\t', skiprows=28)
    read_file_trimmed = read_file.loc[(read_file['IDENTIFIER_TYPE'] == 'SWISS-PROT') &
                                      (read_file['ORGANISM_OFFICIAL_NAME'] == 'Homo sapiens')]
    biogrid_swissprot = read_file_trimmed[['BIOGRID_ID' , 'IDENTIFIER_VALUE']]
    return biogrid_swissprot.to_csv(out_filename, sep='\t', index=False)

if __name__ ==  '__main__':
    get_biogrid_swiss_conv(in_filename='BIOGRID-IDENTIFIERS-3.4.158.tab.txt',
                           out_filename='biogrid_swiss_conversion.csv')