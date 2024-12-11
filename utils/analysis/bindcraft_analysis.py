import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Align import MultipleSeqAlignment

class utils:
    def MSA(ids, seqs):
        assert len(ids) == len(seqs)
        n = len(ids)
        align = MultipleSeqAlignment([])
        for i in range(n):
            align.add_sequence(ids[i], seqs[i])

    def sequence_diversity(seqs):
        pass

    def structure_diversity(structure):
        pass

    def cluster(features):
        pass

    def esm_embed_sequences(seq):
        pass

    def get_biopython_features(seq):
        pass

    def read_fastas(fastas_dir):
        pass

    def write_fastas(ids, seqs):
        pass

class analysis():
    def __init__(self, metrics_filepath, output_dir='', index_column='Design', target_protein_length=[100], target_chains=["A"]):
        self.df = pd.read_csv(metrics_filepath)
        self.target_protein_length = target_protein_length
        self.target_protein_chains = target_protein_chains
        if index_column != '': 
            self.df.index = self.df[index_column]
        self.out_dir = output_dir

        self.boxplot()
        self.corr_matrix()
        self.cov_matrix()
        self.extract_bound_residues()

    def boxplot(self, columns=['PLDDT', 'I_pLDDT', 'PTM', 'I_pTM', 'PAE', 'I_pAE', 'Ss_pLDDT', 'Surface_hydrophobicity', 'Shape Complementarity', 'Pack Stat']):
        self.df.filter(columns).boxplot()
        plt.savefig(self.out_dir+'boxplot.png') 

    def corr_matrix(self, corr_type='spearman'):
        numerics = df.select_dtypes(include=numerics)
        fig, ax = plt.subplot()
        corr_matrix = df.filter(numerics).corr(corr_type)
        sns.heatmap(corr_matrix, xticklabels=numerics.columns, yticklabels=numerics.columns, ax=ax)
        ax.set_title("Correlation")
        plt.savefig(self.out_dir+'corr_matrix.png')

    def cov_matrix(self):
        numerics = df.select_dtypes(include=numerics)

        fig, ax = plt.subplot()
        cov_matrix = df.filter(numerics).cov()
        sns.heatmap(cov_matrix, xticklabels=numerics.columns, yticklabels=numerics.columns, ax=ax)
        ax.set_title("Covariance")
        plt.savefig(self.out_dir+'cov_matrix.png')

    def extract_bound_residues(self):
        labels = []
        for chain in self.target_protein_chains:
            for i in range(self.target_protein_length):
                labels.append(chain+str(i))
        
        n_design = len(self.df['Design'])
        
        binding_matrix = np.zeros(len(labels), n_design)

        for i, design in enumerate(self.df.iterrows()):
            for j, res in enumerate(design['Interface Residue']):
                idx = labels.index(res)
                binding_matrix[j,i] = 1
        
        fig, ax = plt.subplot()
        sns.heatmap(binding_matrix, xticklabels=labels, yticklabels=self.df['Design'], ax=ax)

        plt.savefig(self.out_dir+'bonding_matrix.png')

    def deep_scatter(self, columns):
        n_cols = len(columns)
        fig, ax = plt.subplots(n_cols, n_cols)
        for i, col_i in enumerate(columns):
            for j, col_j in enumerate(columns):
                ax[i,j].scatter(self.df[col_i], self.df[col_j])
        
        for i in range(n_cols):
            ax[i, -1].set_xlabel(columns[i])
            ax[0, i].set_ylabel(columns[i])

        plt.savefig(self.out_dir+'eval_df_scatter.png')

if __name__ == '__main__':
    a = analysis('inputs/IgE_run1_20241101.csv', output_dir='outputs', target_protein_length=[103,103], target_protein_chains=["B","D"])
    


