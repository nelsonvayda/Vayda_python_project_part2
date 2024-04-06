#!/usr/bin/env python
# coding: utf-8

# In[6]:


import os
from Bio import SeqIO
import pandas as pd

# Defines data extraction from genbank files
def extract_data(file_path):
    data = []
    for seq_record in SeqIO.parse(file_path, "genbank"):
        accession = seq_record.id
        num_features = len(seq_record.features)
        family = genus = species = source = 'Unknown'

        # Gets taxonomic information from files
        if 'taxonomy' in seq_record.annotations:
            taxonomy = seq_record.annotations['taxonomy']
            if len(taxonomy) >= 2:
                family = taxonomy[-2]
                genus = taxonomy[-3]
            elif len(taxonomy) == 1:
                genus = taxonomy[-1]

        # Gets source information from files
        for feature in seq_record.features:
            if feature.type == "source":
                if 'organism' in feature.qualifiers:
                    species = feature.qualifiers['organism'][0]
                if 'isolate' in feature.qualifiers or 'clone' in feature.qualifiers:
                    source = feature.qualifiers.get('isolate', feature.qualifiers.get('clone', ['Unknown']))[0]
        
        data.append([accession, family, genus, species, num_features, source])
    return data

def main():
    # Gets the directory where the script is running
    script_dir = '.'
    
    # Lists all genbank files in the directory
    genbank_files = [f for f in os.listdir(script_dir) if f.endswith('.gb') or f.endswith('.genbank')]
    
    all_data = []
    for file_name in genbank_files:
        file_path = os.path.join(script_dir, file_name)
        all_data.extend(extract_data(file_path))
    
    # Compiles data from all genbank files into a DataFrame
    df = pd.DataFrame(all_data, columns=['Accession', 'Family', 'Genus', 'Species', 'Num_Features', 'Source'])
    
    # Defines the output CSV file name
    output_file = os.path.join(script_dir, "genbank_parse.csv")
    df.to_csv(output_file, index=False)
    
    print(f"Data from {len(genbank_files)} files compiled into {output_file}")

if __name__ == "__main__":
    main()


# In[ ]:




