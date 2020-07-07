import pandas as pd
import pickle


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
    return (mhc1, mhc2)

def load_frequency_dfs(path):
    mhc1 = pd.read_pickle(path + 'IEDB_population_frequency2392_normalized.pkl')
    mhc2 = pd.read_pickle(path + 'IEDB_population_frequency_mhc2_275normalized.pkl')
    return (mhc1, mhc2)


### generate dictionary of genotypes associated with each HLA allele
# get set of genotypes associated with given allele
def get_allele_genotypes(allele, frequency_df, genotype_threshold):
    genotypes = set()
    df = frequency_df.xs(allele, axis=1, level=1, drop_level=False)
    df.columns = df.columns.droplevel()
    df = df[df[allele] >= genotype_threshold]
    for genotype in df.index:
        genotypes.add(genotype)
    return genotypes

# iterate through all alleles and generate dictionary
def generate_allele_genotype_dict(frequency_dfs, genotype_threshold):
    (freq_mhc1, freq_mhc2) = frequency_dfs
    alleles_mhc1 = freq_mhc1.columns.get_level_values(1)
    alleles_mhc2 = freq_mhc2.columns.get_level_values(1)
    
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

## main call
# iterate through all peptides & alleles and generate dataframes
def generate_objective_dfs(affinity_dfs, frequency_dfs, allele_genotype_dict, affinity_threshold):
    (affinity_mhc1, affinity_mhc2) = affinity_dfs
    alleles_mhc1 = affinity_mhc1.columns.get_level_values(1)
    alleles_mhc2 = affinity_mhc2.columns.get_level_values(1)
    
    # modify columns to fit desired data
    objective_columns = pd.MultiIndex.from_product([["Genotypes"], get_all_genotypes(frequency_dfs)])
    
    mhc1 = pd.DataFrame(0, index=affinity_mhc1.index, columns=objective_columns)
    mhc2 = pd.DataFrame(0, index=affinity_mhc2.index, columns=objective_columns)
    
    for peptide in affinity_mhc1.index:
        for allele in alleles_mhc1:
            # modify function call to update desired data
            pass
            
    for peptide in affinity_mhc2.index:
        for allele in alleles_mhc2:
            # modify function call to update desired data
            pass            
    
    return (mhc1, mhc2)


### create summary dataframe
def merge_objective_dfs(objective_dfs):
    (objective_mhc1, objective_mhc2) = objective_dfs
    
    objective = pd.concat([objective_mhc1, objective_mhc2])
    
    return objective

def create_summary_df(feature_df, objective_dfs):
    objective = merge_objective_dfs(objective_dfs)

    summary = feature_df.copy(deep=True)
    summary.index = summary.index.rename("Peptide")
    summary.columns = pd.MultiIndex.from_product([["Features"], summary.columns])
    summary = pd.merge(summary, objective, how="inner", on="Peptide")

    return summary


### store summmary dataframe
def store_summary(summary_df, sample, n = 100):    
    # take subset
    if sample:
        summary_df = summary_df.sample(n=n)

    # store as pickle
    summary_df.to_pickle(('summary-%d.pkl' % n) if sample else 'summary.pkl')


### main
def main():
    frequency_dfs = load_frequency_dfs(path)
    allele_genotype_dict = generate_allele_genotype_dict(frequency_dfs, genotype_threshold)
    affinity_dfs = load_affinity_dfs(path)
    objective_dfs = generate_objective_dfs(affinity_dfs, frequency_dfs, allele_genotype_dict, affinity_threshold)
    feature_df = load_feature_df(path)
    summary_df = create_summary_df(feature_df, objective_dfs)
    print(summary_df)
    store_summary(summary_df, sample=False)

if __name__ == '__main__':
    main()