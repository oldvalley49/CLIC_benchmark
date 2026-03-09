import itertools
import os

tissues = ["PBMC", "BMMC-s4d8","BMMC-s1d2","TDBM", "mBrain", "mRetina", "mSkin"]
var_nums = [2000, 4000, 5000, 6000, 8000, 10000, 12000, 14000]
total_nums = [2000]
sub_nums = [1]
indices = range(1,4) 
cor_methods = ["signac-pearson", "maestro-pearson"]
activity_models = ['signac', 'maestro']

combinations = []


tissue_combinations = list(itertools.product(tissues, var_nums, total_nums, sub_nums, cor_methods, indices, activity_models))
combinations.extend(tissue_combinations)

with open("scripts/benchmark/parameters/parameters_paired.txt", "w") as f:
    for combo in combinations:
        f.write(" ".join(map(str, combo)) + "\n")
f.close()

tissues_rna = ["BMMC-s1d2", "BMMC-s4d8", "BMMC-s3d6"]
tissues_atac = ["BMMC-s2d4", "BMMC-s1d2", "BMMC-s2d4"]
cor_methods = ["maestro_pearson", "signac_pearson"]
var_nums =  [2000, 4000, 5000, 6000, 8000, 10000, 12000, 14000]


combinations = []
for i in range(len(tissues_rna)):
    tissue_combinations = list(itertools.product([tissues_rna[i]], [tissues_atac[i]], var_nums, total_nums, cor_methods, indices))
    combinations.extend(tissue_combinations)

with open("scripts/benchmark/parameters/parameters_unpaired.txt", "w") as f:
    for combo in combinations:
        f.write(" ".join(map(str, combo)) + "\n")
f.close()