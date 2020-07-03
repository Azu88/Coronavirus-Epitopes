import pandas as pd
import pickle


## parameters
genotype_threshold = 0.01
affinity_threshold = 0.638


## load relevant data
path = 'C:\\Users\\Nika\\Downloads\\Organized\\Coronavirus-Epitopes\\optivax-master\\'
all_epitope_features = pd.read_pickle(path + 'AllEpitopeFeatures.pkl')
affinity_pred_mhc1 = pd.read_pickle(path + 'dropbox_data\\mean2_mhc1_pred_affinity_pivot.pkl')
affinity_pred_mhc2 = pd.read_pickle(path + 'dropbox_data\\netmhc4.0_mhc2_pred_affinity_pivot.pkl')
genotype_frequency_mhc1 = pd.read_pickle(path + 'IEDB_population_frequency2392_normalized.pkl')
genotype_frequency_mhc2 = pd.read_pickle(path + 'IEDB_population_frequency_mhc2_275normalized.pkl')


## initialize summary dataframe
summary = all_epitope_features.copy(deep=True)
summary.index = summary.index.rename("Peptide")
summary.columns = pd.MultiIndex.from_product([["Features"], summary.columns])
summary = pd.merge(summary, affinity_pred_mhc1, how="left", on="Peptide")
summary = pd.merge(summary, affinity_pred_mhc2, how="left", on="Peptide")
genotypes = pd.Index.union(genotype_frequency_mhc1.index, genotype_frequency_mhc2.index)
genotype_columns = pd.MultiIndex.from_product([["Genotypes_min", "Genotypes_max"], genotypes])
genotype_zeros = pd.DataFrame(0, index=summary.index, columns=genotype_columns)
summary = pd.merge(summary, genotype_zeros, how="left", on="Peptide")


## take subset
summary = summary.sample(n=100)


## create dictionary of genotypes associated with each HLA
# create set of genotypes associated with given allele
def associated_genotypes(allele, frequency_df, genotype_threshold):
    genotype_set = set()
    df = frequency_df.xs(allele, axis=1, level=1, drop_level=False)
    df.columns = df.columns.droplevel()
    df = df[df[allele] >= genotype_threshold]
    for genotype in df.index:
        genotype_set.add(genotype)
    return genotype_set

# iterate through all alleles
genotype_dict = dict()
alleles_mhc1 = genotype_frequency_mhc1.columns.get_level_values(1)
alleles_mhc2 = genotype_frequency_mhc2.columns.get_level_values(1)

for allele in alleles_mhc1:
    genotype_dict[allele] = associated_genotypes(allele, genotype_frequency_mhc1, genotype_threshold)
    
for allele in alleles_mhc2:
    genotype_dict[allele] = associated_genotypes(allele, genotype_frequency_mhc2, genotype_threshold)


## determine min & max affinity for genotypes associated with each peptide
# update min & max affinity for all relevant genotypes given a peptide & allele
def updateMinMax(summary, genotype_dict, peptide, allele, affinity_threshold):
    affinity = summary.loc[peptide, summary.columns.get_level_values(1)==allele].iloc[0]
    if affinity >= affinity_threshold:
        for genotype in genotype_dict[allele]:
            genotype_min = summary.loc[peptide, ("Genotypes_min", genotype)]
            if genotype_min == 0 or affinity < genotype_min:
                summary.at[peptide, ("Genotypes_min", genotype)] = affinity
            genotype_max = summary.loc[peptide, ("Genotypes_max", genotype)]
            if affinity > genotype_max:
                summary.at[peptide, ("Genotypes_max", genotype)] = affinity

# iterate through all peptides and alleles
for peptide in summary.index:
    for allele in alleles_mhc1:
        updateMinMax(summary, genotype_dict, peptide, allele, affinity_threshold)
    for allele in alleles_mhc2:
        updateMinMax(summary, genotype_dict, peptide, allele, affinity_threshold)


## store summary as pickle
summary.to_pickle('summary.pkl')