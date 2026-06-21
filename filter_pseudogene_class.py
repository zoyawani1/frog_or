import pandas as pd 
import numpy as np 
import re 
import os 
from dataclasses import dataclass 


@dataclass 
class FilterPseudoGeneParams():
    
    # data_path = './data/'
    query_file = None
    human_seq = None
    domain_file = 'HumanORTMD.txt' 
    
    # Filter parameters 
    filter_domains = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']
    filter_consecutive_mismatch_threshold = 5

    output_path = './output/'
    save_data = False
    save_fasta = False

class FilterPseudoGene():
    def __init__(self, params: dataclass):
        self.params = params 
        
        # Read in domain reference dataframe 
        self.domain_df = pd.read_csv(self.params.domain_file)
        print(f"Domain file : {self.params.domain_file}")
        
        # Read in fasta file to be queried into dataframe
        self.query_df = self.read_fasta(self.params.query_file)
        print(f"Query file : {self.params.query_file}")

        # Extract the Human sequence as reference from query_df 
        self.human_seq = self.query_df[self.query_df['name'].str.contains('Human')].sequence[0]
        self.domain_df = self.get_domain_info(self.domain_df, 
                                              human_seq = self.human_seq)
    
    def filter_genes(self):
       # os.makedirs(self.params.output_path, exist_ok=True)
        # Initialize mismatch as False until mismatches are found 
        self.query_df['mismatch'] = False
        # Iterate through rows (sequence entries) and determine for mismatches against Human_seq 
        for index, row in self.query_df.iterrows():
            for domain in self.params.filter_domains:
                start = int(self.domain_df[self.domain_df['domain'] == domain].start)
                end = int(self.domain_df[self.domain_df['domain'] == domain].end)

                # Searching for mismatches between seq 
                mismatch = self.seq_mismatch(self.human_seq[start:end],
                                             row.sequence[start:end], 
                                             consecutive_mismatch_threshold = self.params.filter_consecutive_mismatch_threshold)

                if mismatch:
                    self.query_df.loc[self.query_df['name'] == row['name'], 'mismatch'] = mismatch
                    break
                
                # If there are < 5 aa after the TM7 domain then mismatch also True
                if domain == 'TM7':
                    c_term_len = len(row.sequence[end+1:len(row.sequence)].replace('-',''))
                    if c_term_len <= self.params.filter_consecutive_mismatch_threshold:
                        self.query_df.loc[self.query_df['name'] == row['name'], 'mismatch'] = True
        
        
        if self.params.save_data:
            # Save entire query_df to csv 
            file_save_path = os.path.join(self.params.output_path, os.path.basename(self.params.query_file).split('.fasta')[0])+'_filtered.csv'
            print(f"Saving query df as csv: {file_save_path}")
            self.query_df.to_csv(file_save_path, 
                                 index=0)
        if self.params.save_fasta:
            file_save_path = os.path.join(self.params.output_path, os.path.basename(self.params.query_file).split('.fasta')[0])+'_filtered.fasta'
            print(f"Saving filtered query fasta: {file_save_path}")
            # Saving non-mismatched sequence to fasta
            with open(file_save_path, 'w') as file:
                for _, row in self.query_df.iterrows():
                    if not row['mismatch']:
                        file.write(f">{row['name']}\n")
                        file.write(f"{row['sequence']}\n")
                        
    
    
    def read_fasta(self, filename):
        data = []
        with open(filename, 'r') as file:
            lines = file.readlines()
            sequence_name = None
            sequence = ''
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    if sequence_name is not None:
                        data.append((sequence_name, sequence))
                        sequence = ''
                    sequence_name = line[1:]
                else:
                    sequence += line
            if sequence_name is not None:
                data.append((sequence_name, sequence))
        data = pd.DataFrame(data, columns = ['name', 'sequence']) 
        return data

    def get_domain_info(self, 
                        domain_reference: pd.DataFrame,
                        human_seq:str 
                        ):

        for q_motif in domain_reference.motif:
            
            pattern = '-*'+('-*'.join(list(q_motif)))
            matches = re.finditer(pattern, human_seq)
            
            for match in matches:
                start_pos = match.start()
                end_pos = match.end() - 1  # Adjust for inclusive end position    
            
                
            domain_reference.loc[domain_reference['motif'] == q_motif,['start']] = start_pos
            domain_reference.loc[domain_reference['motif'] == q_motif,['end']] = end_pos
        
        
        domain_reference[['start', 'end']] = domain_reference[['start', 'end']].astype(int)
        
        return domain_reference

        
    def seq_mismatch(self, 
                     ref_seq: str, 
                     q_seq:str, 
                     consecutive_mismatch_threshold = 5
                     ):
        
        consecutive_mismatch_count = 0
        mismatch = False
        for ref_char, q_char in zip(ref_seq, q_seq):
            if (ref_char == '-') and (q_char == '-'):
                continue

            if (ref_char != '-' and q_char == '-') or (ref_char == '-' and q_char != '-'):
                consecutive_mismatch_count += 1
                if (consecutive_mismatch_count >= consecutive_mismatch_threshold) | (q_char == '*'):
                    mismatch = True    
                    break
            else:
                consecutive_mismatch_count = 0

        return mismatch
            
