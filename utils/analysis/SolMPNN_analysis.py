import re
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqIO import FastaIO
from Bio.Cluster import kcluster#

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib.pyplot as plt
import seaborn as sns

import torch
import esm

import argparse

def embed_sequence(sequences):
        
    """Load ESM-2 model"""
    print(len(sequences), sequences)
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    batch_labels, batch_strs, batch_tokens = batch_converter(sequences)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]
    print(token_representations)
    # Sequence representation
    # sequence_representations = []
    # for i, tokens_len in enumerate(batch_lens):
    #     sequence_representations.append(token_representations[i, 1 : tokens_len - 1])
    
    return token_representations

def filter_sequences(input_file, output_file): #TODO: rmeove file writing and making it return and add a file writing with it.
    print('filtering')
    print(input_file, output_file)
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        entry_name = ""
        sequence = ""
        write_entry = False

        for line in infile:
            if line.startswith('>'):
                # Write the previous entry if it's valid
                if entry_name and sequence and not has_repeated_chars(sequence):
                    outfile.write(entry_name + '\n')
                    outfile.write(sequence + '\n')
                # Start a new entry
                entry_name = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        
        # Write the last entry if it's valid
        if entry_name and sequence and not has_repeated_chars(sequence):
            print('cringe')
            outfile.write(entry_name + '\n')
            outfile.write(sequence + '\n')

# Function to extract values from the sequence description
def extract_values(description):
    try:
        parts = description.split(', ')
        seq_rec = float(parts[6].split('=')[1])
        overall_confidence = float(parts[4].split('=')[1])
        ligand_confidence = float(parts[5].split('=')[1])
        return seq_rec, overall_confidence, ligand_confidence
    except (IndexError, ValueError) as e:
        print(f"Error parsing description: {description}")
        print(f"Error: {e}")
        return None, None, None

def has_repeated_chars(sequence, max_repeats=3):
    pattern = re.compile(r"(.)\1{" + str(max_repeats) + ",}")
    return pattern.search(sequence) is not None

def score_sequences(fasta_file):
    # List to store results
    results = []
    seq_records = []  # List to store sequence records

    # Parse the FASTA file and compute protein parameters for each sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.description.split(', ')[1].split('=')[1]  # Extract the id
        seq_rec, overall_confidence, ligand_confidence = extract_values(record.description)
                
        if seq_rec is None or overall_confidence is None or ligand_confidence is None:
            continue  # Skip this record if there was an error

        seq = str(record.seq)
        analysis = ProteinAnalysis(seq)
        
        # Compute various protein parameters
        pi = analysis.isoelectric_point()
        molecular_weight = analysis.molecular_weight()
        aromaticity = analysis.aromaticity()
        instability_index = analysis.instability_index()
        gravy = analysis.gravy()
        secondary_structure_fraction = analysis.secondary_structure_fraction()
        amino_acid_percent = analysis.get_amino_acids_percent()
        amino_acid_count = analysis.count_amino_acids()
        extinction_coefficient = analysis.molar_extinction_coefficient()

        # Append results to the list
        results.append({
            "ID": seq_id,
            "Seq Rec": seq_rec,
            "Overall Confidence": overall_confidence,
            "Ligand Confidence": ligand_confidence,
            "Isoelectric Point": pi,
            "Molecular Weight": molecular_weight,
            "Aromaticity": aromaticity,
            "Instability Index": instability_index,
            "GRAVY": gravy,
            "Helix Fraction": secondary_structure_fraction[0],
            "Turn Fraction": secondary_structure_fraction[1],
            "Sheet Fraction": secondary_structure_fraction[2],
            "Extinction Coefficient (reduced)": extinction_coefficient[0],
            "Extinction Coefficient (oxidized)": extinction_coefficient[1],
            "seq" : seq
        })
        
        # Store the sequence record for later use
        seq_records.append(record)

    # Create a DataFrame
    df = pd.DataFrame(results)
    print("DataFrame columns:", df.columns)  # Debug statement to check columns

    # Find min and max for normalization
    overall_conf_min = df['Overall Confidence'].min()
    overall_conf_max = df['Overall Confidence'].max()
    ligand_conf_min = df['Ligand Confidence'].min()
    ligand_conf_max = df['Ligand Confidence'].max()

    # Function to normalize values
    def normalize(value, min_val, max_val):
        return (value - min_val) / (max_val - min_val)

    # # Normalize the "Overall Confidence" and "Ligand Confidence" columns
    # df['Normalized Overall Confidence'] = df['Overall Confidence'].apply(normalize, args=(overall_conf_min, overall_conf_max))
    # df['Normalized Ligand Confidence'] = df['Ligand Confidence'].apply(normalize, args=(ligand_conf_min, ligand_conf_max))

    # # Calculate the average of the normalized values
    # df['Average Normalized Confidence'] = df[['Normalized Overall Confidence', 'Normalized Ligand Confidence']].mean(axis=1)

    # Calculate the average for each numeric column and add it as a new row
    average_row = {col: df[col].mean() if df[col].dtype in ['float64', 'int64'] else 'Average' for col in df.columns}

    # Convert the average_row dictionary into a DataFrame
    

    output_file = fasta_file.replace('.fa', '_conf_protparam_normalized.csv')

    # Save to CSV
    df.to_csv(output_file, index=False)
    return df

def top_n(df, n, output_file='topn.fa'):
        
    # Get the top 10 sequences based on "Average Normalized Confidence"
    top_10_df = df.nlargest(n, 'Average Normalized Confidence')

    # Create a new FASTA file with the top 10 sequences
    top_10_fasta_file = fasta_file.replace('.fa', '_top_10.fa')

    with open(top_10_fasta_file, 'w') as output_handle:
        fasta_writer = FastaIO.FastaWriter(output_handle, wrap=None)  # Set wrap=None to avoid line breaks
        fasta_writer.write_header()
        for seq_id in top_10_df['ID']:
            for record in seq_records:
                if record.description.split(', ')[1].split('=')[1] == seq_id:
                    fasta_writer.write_record(record)
                    break
        fasta_writer.write_footer()

    # Print completion message with the file path
    print("Analysis complete. Results saved to " + output_file)
    print("Top 10 sequences saved to " + top_10_fasta_file)

class SolMPNN_run:
    def __init__(self, fasta_filepath, process=True, top_n=5, out_dir='', chain_index=-1):
        self.fasta_filepath = fasta_filepath
        self.out_dir = out_dir
        self.chain_index = chain_index

        self.load_fasta_file(process)

        self.df = score_sequences(self.fasta_filepath_processed)

        self.top_n_OC = self.df.sort_values("Overall Confidence", ascending=False).iloc[:top_n]
        self.top_n_gravy = self.df.sort_values("GRAVY", ascending=False).iloc[:top_n]
        print()


        self.df_analysis()
        # print(self.df.iloc[0], self.df.iloc[-1], self.df, len(self.df.index))
        self.protein_dataset = list(zip(self.df['ID'], self.df['seq']))
        self.df['embedded_sequence'] = embed_sequence(self.protein_dataset)

        self.cluster_analysis()
    
    def cluster_analysis(self, n_components=2, verbose=0, perplexity=40, n_iter=300):
        pca_instance = PCA(n_components=n_components)
        print(self.df['embedded_sequence'].values.shape)
        pca_results = pca_instance.fit_transform(self.df['embedded_sequence'].values)

        tsne = TSNE(n_components=n_components)
        tSNE_results = tsne.fit_transform(pca_results)

        plt.figure(figsize=(16,10))
        sns.scatterplot(
            x=tSNE_results[:,0], y=tSNE_results[:,1],
            color=self.df['Normalized Overall Confidence'],
            legend="full",
            alpha=0.3
        )
        plt.savefig(self.out_dir+'/ESM_tSNE_sequence.png')
        
    def df_analysis(self, corr_type='spearman'):
        numerics = self.df.select_dtypes(include='number')

        corr_matrix = numerics.corr(corr_type)
        sns.heatmap(corr_matrix)
        plt.subplots_adjust(left=0.3, bottom=0.3)
        plt.savefig(self.out_dir+'Corr_matrix.png')

        cov_matrix = numerics.cov()
        sns.heatmap(cov_matrix)
        plt.subplots_adjust(left=0.3, bottom=0.3)

        plt.savefig(self.out_dir+'Cov_matrix.png')

        bp_df = numerics.drop(columns=['GRAVY', 'Isoelectric Point', 'Instability Index', 'Molecular Weight', 'Extinction Coefficient (reduced)', 'Extinction Coefficient (oxidized)'])
        sns.boxplot(bp_df, xlabel=bp_df.columns)
        plt.savefig(self.out_dir+'boxplot.png')

    def load_fasta_file(self, process):
        with open(self.fasta_filepath, 'r') as f:
            text = f.read()
        
        if process:
            self.fasta_text = self.process_fasta_file(text)
        else:
            self.fasta_text = text
            self.fasta_filepath_processed = self.fasta_filepath

    def process_fasta_file(self, text):
        # print(text.split('\n')[2:])
        self.input_fasta_text = text.split('\n')[:2]
        design_fasta_text = ''
        for line in text.split('\n')[2:]:
            if not line.startswith('>'):
                design_fasta_text += line.split(':')[self.chain_index] + '\n'
            else:
                design_fasta_text += line + '\n'

        self.fasta_filepath_processed = self.fasta_filepath[:-3]+'_processed.fa'
        with open(self.fasta_filepath, 'w+') as f:
            f.write(design_fasta_text)
        
        filter_sequences(self.fasta_filepath, self.fasta_filepath_processed)
        
        with open(self.fasta_filepath_processed, 'r') as f:
            print(f.read())
        

        return design_fasta_text



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to run BindCraft binder design.')
    parser.add_argument('--fasta_filepath')
    parser.add_argument('--output_dir', default='/user/home/vi21227/code/vi21227/code/ProteinDesign/outputs/solMPNN/')
    parser.add_argument('--chain_idx', default=0)
    args = parser.parse_args()
    print('inps', args.fasta_filepath[:-4]+'.fa', args.output_dir)
    solmpnn = SolMPNN_run(args.fasta_filepath[:-4]+'.fa', out_dir=args.output_dir, chain_index=int(args.chain_idx))

