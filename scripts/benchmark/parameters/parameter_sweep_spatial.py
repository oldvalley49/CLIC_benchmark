import itertools
import os


tissues = ["mBrain", "mEmbryo2", "Melanoma"]
var_nums = [2000, 4000, 5000, 6000, 8000]
total_nums = [2000]
sub_nums = [1]
indices = range(1, 4)  


combinations = []

for tissue in tissues:
    if tissue.startswith('m'):
        cor_methods = ["corr"]
    else:
        cor_methods = ["corr"]
    
    tissue_combinations = list(itertools.product([tissue], var_nums, total_nums, sub_nums, cor_methods, indices))
    combinations.extend(tissue_combinations)

# save combinations to a file
with open("scripts/benchmark/parameters/parameters_paired_spatial.txt", "w") as f:
    for combo in combinations:
        f.write(" ".join(map(str, combo)) + "\n")
f.close()
