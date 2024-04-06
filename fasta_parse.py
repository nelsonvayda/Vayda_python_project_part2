#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
from Bio import SeqIO
import pandas as pd

# Defines function to extract data from fasta file
def extract_data(file_path):
    data = []
    for seq_record in SeqIO.parse(file_path, "fasta"):
        fasta_id = seq_record.id
        first_10_aa = seq_record.seq[:10]
        length_of_protein = len(seq_record.seq)
        total_c = seq_record.seq.count('C')
        
        data.append([fasta_id, str(first_10_aa), length_of_protein, total_c])
    return data


def main():
    # Gets directory where script is running
    script_dir = '.'

    fasta_files = [f for f in os.listdir(script_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    
    all_data = []
    for file_name in fasta_files:
        file_path = os.path.join(script_dir, file_name)
        all_data.extend(extract_data(file_path))
    
    # Compiles data from all files into a DataFrame
    df = pd.DataFrame(all_data, columns=['ID', 'First_10_AA', 'Length', 'Number_Cs'])
    
    output_file = os.path.join(script_dir, "protein_info.csv")
    df.to_csv(output_file, index=False)
    
    print(f"Data from {len(fasta_files)} files compiled into {output_file}")

if __name__ == "__main__":
    main()


# In[ ]:




