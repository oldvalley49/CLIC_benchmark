import itertools
import os

# Define parameters
tissues = ["PBMC", "BMMC-s4d8","BMMC-s1d2","TDBM", "mBrain", "mSkin", "mRetina"]
var_nums = [16000, 18000]
total_nums = [2000]
sub_nums = [1]
indices = range(1,4)  # 3 iterations per condition
cor_methods = ["corr"]

# Generate all combinations
combinations = []

    
# Generate combinations for the current tissue
tissue_combinations = list(itertools.product(tissues, var_nums, total_nums, sub_nums, cor_methods, indices))
combinations.extend(tissue_combinations)

# Save combinations to a file
with open("scripts/benchmark/parameters/parameters_paired.txt", "w") as f:
    for combo in combinations:
        f.write(" ".join(map(str, combo)) + "\n")
f.close()

tissues_rna = ["BMMC-s1d2", "BMMC-s4d8", "BMMC-s3d6"]
tissues_atac = ["BMMC-s2d4", "BMMC-s1d2", "BMMC-s2d4"]
cor_methods = ["corr"]
var_nums = [16000, 18000]


combinations = []
for i in range(len(tissues_rna)):
    tissue_combinations = list(itertools.product([tissues_rna[i]], [tissues_atac[i]], var_nums, total_nums, cor_methods, indices))
    combinations.extend(tissue_combinations)

with open("scripts/benchmark/parameters/parameters_unpaired.txt", "w") as f:
    for combo in combinations:
        f.write(" ".join(map(str, combo)) + "\n")
f.close()