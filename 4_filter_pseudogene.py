
import numpy as np 
import pandas as pd 
import sys
import os
import argparse

from filter_pseudogene_class import FilterPseudoGeneParams, FilterPseudoGene


parser = argparse.ArgumentParser()
parser.add_argument('--save_data', 
                    help = "boolean to save data", 
                    default=False)
parser.add_argument('--save_fasta', 
                    help = "boolean to save data", 
                    default=True)
parser.add_argument('--output_path', 
                    help = "file path to output data", 
                    default='./')
parser.add_argument('--query_file_path', 
                    help = "file path of .fasta to be queried ", 
                    required=True)
parser.add_argument('--domain_file_path', 
                    help = "file path of domain info", 
                    default = 'HumanORTMD.txt')
parser.add_argument('-get_tm1_pos',  # Add the -get_domain flag as an optional argument
                    help="Flag to indicate getting the first amino acid of aligned TM1",
                    action="store_true")  # This makes it a boolean flag
args = parser.parse_args()


# instantiate parameters
filter_params = FilterPseudoGeneParams()   
filter_params.save_data = args.save_data 
filter_params.save_fasta = args.save_fasta
filter_params.output_path = args.output_path
filter_params.query_file = args.query_file_path
filter_params.domain_file = args.domain_file_path


# Check if the -get_domain flag is provided
if args.get_tm1_pos:
    print("===== Getting domain =====\n")
    filter_OR = FilterPseudoGene(filter_params)
    tm1_position = int(filter_OR.domain_df[filter_OR.domain_df['domain'] == 'TM1'].start.values)+1
    print(f'First amino acid at TM1 is : {tm1_position}')
    sys.exit()  # Stop the script after printing

# instantiate class and run filter_genes 
filter_OR = FilterPseudoGene(filter_params)
filter_OR.filter_genes()
