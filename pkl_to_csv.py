import pandas as pd
import pickle

path = 'C:\\Users\\Nika\\Downloads\\Organized\\Coronavirus-Epitopes\\optivax-master\\'
all_epitope_features = pd.read_pickle(path + 'AllEpitopeFeatures.pkl')
affinity_pred_mhc1 = pd.read_pickle(path + 'dropbox_data\\mean2_mhc1_pred_affinity_pivot.pkl')
affinity_pred_mhc2 = pd.read_pickle(path + 'dropbox_data\\netmhc4.0_mhc2_pred_affinity_pivot.pkl')
genotype_frequency_mhc1 = pd.read_pickle(path + 'IEDB_population_frequency2392_normalized.pkl')
genotype_frequency_mhc2 = pd.read_pickle(path + 'IEDB_population_frequency_mhc2_275normalized.pkl')

all_epitope_features.to_csv('AllEpitopeFeatures.csv')
affinity_pred_mhc1.to_csv('mean2_mhc1_pred_affinity_pivot.csv')
affinity_pred_mhc2.to_csv('netmhc4.0_mhc2_pred_affinity_pivot.csv')
genotype_frequency_mhc1.to_csv('IEDB_population_frequency2392_normalized.csv')
genotype_frequency_mhc2.to_csv('IEDB_population_frequency_mhc2_275normalized.csv')