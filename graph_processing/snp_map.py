"""
Generate a dictonary to map from SNP id to its postion (a tuple containing node_id, chromosome, chunk_start, chunk_end).
For example: {'rs2489000': (1, 1000000, 1005000, 1)}
"""
import argparse
import pandas as pd
from os import listdir
from os.path import isfile, join
import json

def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-snps_dir',
                        default='../../final_results_5kb_to_Dl/snps_splitted_chr_final/',
                        required=False,
                        help='directory of the 23 SNPs files.')

    parser.add_argument('-snps_map_dir',
                        default='../../snp_map/snp_map.json',
                        required=False,
                        help='directory of the output snp mapping file.')  
                 
    return parser.parse_args()


def load_snps(snps_dir):
    """
    Load SNPs data of the 23 chromosomes into the graph
    """
    f_list = [join(snps_dir, f) for f in listdir(snps_dir) if isfile(join(snps_dir, f))] # get the file list
    f_list = sorted(f_list, key=lambda x: int(x.split('_')[-1])) # sort the list according to suffix
    assert len(f_list) == 23
    
    snp_map = {}
    for i in range(len(f_list)):
        print('extrating info from chromosome {}...'.format(i))
        df_chr = pd.read_csv(f_list[i], delim_whitespace=True)
        snp_map.update(extract_map_from_df(df_chr))
    return snp_map

def extract_map_from_df(df):
    """
    Extract SNP mappings from a dataframe.
    return: a dictionary containing the mapping
    """
    df = df[['Interactor', 'Nodes', 'SNPs_ID']]
    df_new = pd.DataFrame()
    df_new['Nodes'] = df['Nodes']
    df_new['Interactor'] = df['Interactor'].apply(lambda x: [int(e) for e in x.split('-')])
    df_new['SNPs_ID'] = df['SNPs_ID'].apply(lambda x: x.split('__') if type(x)==str else 'nan')
    
    snp_map = {} # dictionary to store the location of the SNPs
    for index, row in df_new.iterrows():
        Interactor = row['Interactor']
        SNPs_ID = row['SNPs_ID']
        Interactor.append(row['Nodes'])
        snp_map.update(dict.fromkeys(SNPs_ID, tuple(Interactor)))
    return snp_map    


if __name__=="__main__":
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 100)

    args = get_args()
    snps_dir = args.snps_dir
    snps_map_dir = args.snps_map_dir
    snp_map = load_snps(snps_dir)

    with open(snps_map_dir, 'w') as outfile:
        json.dump(snp_map, outfile)

