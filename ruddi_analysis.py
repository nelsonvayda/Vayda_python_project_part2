#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
from Bio import SeqIO
import pandas as pd

# Defines function to calculate gc content
def calculate_gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

# Defines function to get ATG codons
def count_atg_codons(seq):
    return seq.count("ATG")

# Defines function to get reverse strand
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(seq)])

def extract_data(file_path):
    data = []
    for seq_record in SeqIO.parse(file_path, "fasta"):
        genome_length = len(seq_record.seq)
        gc_content = calculate_gc_content(seq_record.seq)
        forward_atg_count = count_atg_codons(seq_record.seq)
        reverse_atg_count = count_atg_codons(reverse_complement(seq_record.seq))
        
        data.append([genome_length, gc_content, forward_atg_count, reverse_atg_count])
    return data

def main():
    # Gets the directory where the script is running
    script_dir = '.'
    
    fna_files = [f for f in os.listdir(script_dir) if f.endswith('.fna')]
    
    all_data = []
    for file_name in fna_files:
        file_path = os.path.join(script_dir, file_name)
        all_data.extend(extract_data(file_path))

    # Compiles data from all genbank files into a DataFrame
    df = pd.DataFrame(all_data, columns=['Genome_Length', 'GC_Content', 'Forward_ATG_Count', 'Reverse_ATG_Count'])

    # Defines the output CSV file name
    output_file = os.path.join(script_dir, "ruddi.csv")
    df.to_csv(output_file, index=False)
    
    print(f"Data from {len(fna_files)} files compiled into {output_file}")


if __name__ == "__main__":
    main()


# In[ ]:




