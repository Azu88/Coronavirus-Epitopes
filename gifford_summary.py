import pandas as pd
import pickle

import time # FOR TESTING


### notes
# alleles should be consistent between affinity & frequency data


### parameters
genotype_threshold = 0.01
affinity_threshold = 0.638

path = 'C:\\Users\\Nika\\Downloads\\Organized\\Coronavirus-Epitopes\\optivax-master\\'


### load relevant dataframes
def load_feature_df(path):
    return pd.read_pickle(path + 'AllEpitopeFeatures.pkl')

def load_affinity_dfs(path):
    mhc1 = pd.read_pickle(path + 'dropbox_data\\mean2_mhc1_pred_affinity_pivot.pkl')
    mhc2 = pd.read_pickle(path + 'dropbox_data\\netmhc4.0_mhc2_pred_affinity_pivot.pkl')
    # convert all values from float64 to float32
    mhc1 = mhc1.apply(pd.to_numeric, downcast="float")
    mhc2 = mhc2.apply(pd.to_numeric, downcast="float")
    return (mhc1, mhc2)

def load_frequency_dfs(path):
    mhc1 = pd.read_pickle(path + 'IEDB_population_frequency2392_normalized.pkl')
    mhc2 = pd.read_pickle(path + 'IEDB_population_frequency_mhc2_275normalized.pkl')
    # convert all values from float64 to float32
    mhc1 = mhc1.apply(pd.to_numeric, downcast="float")
    mhc2 = mhc2.apply(pd.to_numeric, downcast="float")
    return (mhc1, mhc2)


### generate dictionary of genotypes associated with each HLA allele
# get list of genotypes associated with given allele
def get_allele_genotypes(allele, frequency_df, genotype_threshold):
    genotypes = []
    frequency_df = frequency_df[frequency_df[allele] >= genotype_threshold]
    for genotype in frequency_df.index:
        genotypes.append(genotype)
    return genotypes

# iterate through all alleles and generate allele-genotype dictionary
def generate_allele_genotype_dict(frequency_dfs, genotype_threshold):
    (freq_mhc1, freq_mhc2) = frequency_dfs
    freq_mhc1.columns = freq_mhc1.columns.droplevel()
    freq_mhc2.columns = freq_mhc2.columns.droplevel()
    freq_mhc1 = freq_mhc1.drop(columns=["unknown"])
    freq_mhc2 = freq_mhc2.drop(columns=["unknown"])
    alleles_mhc1 = freq_mhc1.columns
    alleles_mhc2 = freq_mhc2.columns
    
    allele_genotype_dict = dict()

    for allele in alleles_mhc1:
        allele_genotype_dict[allele] = get_allele_genotypes(allele, freq_mhc1, genotype_threshold)
        
    for allele in alleles_mhc2:
        allele_genotype_dict[allele] = get_allele_genotypes(allele, freq_mhc2, genotype_threshold)
        
    return allele_genotype_dict
    

### generate objective dataframes
# get set of all genotypes as index object
def get_all_genotypes(frequency_dfs):
    (freq_mhc1, freq_mhc2) = frequency_dfs
    all_genotypes = pd.Index.union(freq_mhc1.index, freq_mhc2.index)
    return all_genotypes

## options for data-generating functions
# # update min & max affinity for all relevant genotypes given a peptide & allele
# def updateMinMax(summary, genotype_dict, peptide, allele, affinity_threshold):
#     affinity = summary.loc[peptide, summary.columns.get_level_values(1)==allele].iloc[0]
#     if affinity >= affinity_threshold:
#         print(genotype_dict[allele])
#         for genotype in genotype_dict[allele]:
#             genotype_min = summary.loc[peptide, ("Genotypes_min", genotype)]
#             if genotype_min == 0 or affinity < genotype_min:
#                 summary.at[peptide, ("Genotypes_min", genotype)] = affinity
#             genotype_max = summary.loc[peptide, ("Genotypes_max", genotype)]
#             if affinity > genotype_max:
#                 summary.at[peptide, ("Genotypes_max", genotype)] = affinity

# update peptide-genotype values to 1 if genotype is covered by peptide, 0 if not
def update_genotype_cover(allele, peptide, affinity_df, objective_df, allele_genotype_dict, affinity_threshold):
    affinity = affinity_df.at[peptide, allele]
    if affinity >= affinity_threshold:
        for genotype in allele_genotype_dict[allele]:
            objective_df.at[peptide, ("Genotypes", genotype)] = 1

## main call
# iterate through all peptides & alleles and generate objective dataframes
def generate_objective_dfs(affinity_dfs, frequency_dfs, allele_genotype_dict, affinity_threshold):
    (affinity_mhc1, affinity_mhc2) = affinity_dfs
    affinity_mhc1.columns = affinity_mhc1.columns.droplevel()
    affinity_mhc2.columns = affinity_mhc2.columns.droplevel()
    affinity_mhc1 = affinity_mhc1.drop(columns=["unknown"])
    affinity_mhc2 = affinity_mhc2.drop(columns=["unknown"])
    
    # affinity_mhc1 = affinity_mhc1.sample(n=100) # FOR TESTING
    # affinity_mhc2 = affinity_mhc2.sample(n=100) # FOR TESTING
    
    alleles_mhc1 = affinity_mhc1.columns
    alleles_mhc2 = affinity_mhc2.columns
    
    # modify columns to fit desired data
    objective_columns = pd.MultiIndex.from_product([["Genotypes"], get_all_genotypes(frequency_dfs)])
    
    mhc1 = pd.DataFrame(0, index=affinity_mhc1.index, columns=objective_columns), dtype="int32")
    mhc2 = pd.DataFrame(0, index=affinity_mhc2.index, columns=objective_columns), dtype="int32")
    
    peptides = 1 # FOR TESTING
    t = time.process_time() # FOR TESTING
    
    for peptide in affinity_mhc1.index:
        for allele in alleles_mhc1:
            # modify function call to update desired data
            update_genotype_cover(allele, peptide, affinity_mhc1, mhc1, allele_genotype_dict, affinity_threshold)
            
        if peptides % 100 == 0: # FOR TESTING
            print(peptides) # FOR TESTING
            print(time.process_time() - t) # FOR TESTING
        peptides += 1 # FOR TESTING
            
    for peptide in affinity_mhc2.index:
        for allele in alleles_mhc2:
            # modify function call to update desired data
            update_genotype_cover(allele, peptide, affinity_mhc2, mhc2, allele_genotype_dict, affinity_threshold)            
            
        if peptides % 100 == 0: # FOR TESTING
            print(peptides) # FOR TESTING
            print(time.process_time() - t) # FOR TESTING
        peptides += 1 # FOR TESTING
    
    return (mhc1, mhc2)


### create summary dataframe
# merge objective dataframes for mhc1 and mhc2
def merge_objective_dfs(objective_dfs):
    (objective_mhc1, objective_mhc2) = objective_dfs
    
    objective = pd.concat([objective_mhc1, objective_mhc2])
    
    return objective

# merge feature and objective dataframes
def create_summary_df(feature_df, objective_dfs):
    objective = merge_objective_dfs(objective_dfs)

    summary = feature_df.copy(deep=True)
    summary.index = summary.index.rename("Peptide")
    summary.columns = pd.MultiIndex.from_product([["Features"], summary.columns])
    summary = pd.merge(summary, objective, how="inner", on="Peptide")

    return summary


### store summmary dataframe
# write summary dataframe to a file
def store_summary(path, summary_df):    
    # store as pickle
    summary_df.to_pickle(path + 'summary.pkl')


### main
def main():
    frequency_dfs = load_frequency_dfs(path)
    allele_genotype_dict = generate_allele_genotype_dict(frequency_dfs, genotype_threshold)
    affinity_dfs = load_affinity_dfs(path)
    objective_dfs = generate_objective_dfs(affinity_dfs, frequency_dfs, allele_genotype_dict, affinity_threshold)
    feature_df = load_feature_df(path)
    summary_df = create_summary_df(feature_df, objective_dfs)
    print(summary_df)
    store_summary(path, summary_df)

if __name__ == '__main__':
    main()