"""
Generate a dictonary to map from SNP id to its postion (a tuple containing chromosome, chunk_start, chunk_end and node_id).
For example: {'rs2489000': (1, 1000000, 1005000, 1)}
"""
import argparse
import pandas as pd
from os import listdir
from os.path import isfile, join

def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-snps_dir',
                        default='../../final_results_5kb_to_Dl/snps_splitted_chr_final/',
                        required=False,
                        help='directory of output edge file.')

    parser.add_argument('-snps_map_dir',
                        default='../../final_results_5kb_to_Dl/snps_map.json/',
                        required=False,
                        help='directory of output edge file.')  
                 
    return parser.parse_args()


def load_snps(snps_dir):
    """
    Load SNPs data of the 23 chromosomes into the graph
    """
    f_list = [join(snps_dir, f) for f in listdir(snps_dir) if isfile(join(snps_dir, f))] # get the file list
    f_list = sorted(f_list, key=lambda x: int(x.split('_')[-1])) # sort the list according to suffix
    assert len(f_list) == 23
    #print(f_list)
    #df_list = [pd.read_csv(snps_file, delim_whitespace=True) for snps_file in f_list] # list of dataframes
    #df_chr_0 = pd.read_csv(f_list[0], delim_whitespace=True)
    df_chr_0 = pd.read_csv(f_list[0], delim_whitespace=True)

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 100)
    dict_0 = extract_map_from_df(df_chr_0)


def extract_map_from_df(df):
    """
    Extract SNP mappings from a dataframe.
    return: a dictionary containing the mapping
    """
    df = df[['Interactor', 'Nodes', 'SNPs_ID']]
    df['Interactor'] = df['Interactor'].apply(lambda x: [int(e) for e in x.split('-')])
    df['SNPs_ID'] = df['SNPs_ID'].apply(lambda x: x.split('__') if type(x)==str else 'nan')
    print(df)
    
    snp_map = {}
    return snp_map    


if __name__=="__main__":
    args = get_args()
    snps_dir = args.snps_dir
    load_snps(snps_dir)
