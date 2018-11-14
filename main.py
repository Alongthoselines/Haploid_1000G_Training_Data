import numpy as np
import pandas as pd

all_chromosomes=[]
for i in range(1,23):
    all_chromosomes.append(pd.read_csv('./data/chr{}.phased.vcf'.format(i), delim_whitespace=True))
data=pd.concat(all_chromosomes)

# we get the IDs of individuals from previous training scheme with 15 subpopulations
individus=pd.read_csv('./data/labels_of_individuals_training.txt', delim_whitespace=True)['FID'].tolist()

#we keep only important columns 
keep_columns=np.concatenate((data.columns[:9], individus))

#checking that we really got rid of useless data
print(data.shape)
print(data[keep_columns].shape)


#to create new name for half individuals
def new_name(old_name, modif='New'):
    return old_name+modif



haplotypes=data[data.columns[:9]]
haplotypes['SNP_ID']='chr'+haplotypes['CHROM'].map(str)+'_'+'pos'+haplotypes['POS'].map(str)+'_'+haplotypes[REF]
for ind in individus:
    print(ind)
    haplotypes['{}'.format(new_name(ind,'a'))]=data[ind].str.split('|', expand=True)[0]
    haplotypes['{}'.format(new_name(ind,'b'))]=data[ind].str.split('|', expand=True)[1]

# we can find useless columns like this: haplotypes.columns[:9]
# and it is also important to save the final file with as many information columns as the input format requires
#even if the values in these columns are not important

unenssential_columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

# we don't keep some columns
haplotypes=haplotypes.drop(unenssential_columns, axis=1).T

#but because of required format we add some
haplotype.insert(1, 'IID', 0, allow_duplicates=False)
haplotype.insert(2, 'PAT', 0, allow_duplicates=False)
haplotype.insert(3, 'MAT', 0, allow_duplicates=False)
haplotype.insert(4, 'SEX', 0, allow_duplicates=False)
haplotype.insert(5, 'PHENOtype', -9, allow_duplicates=False)



haplotypes.to_csv('haplotypes_1000G.txt', sep='\t', header=True, index=False)
