import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
import numpy as np

tissues = [i for i in os.listdir('output/paired/rna_atac') if not i.startswith('.') and i != 'mRetina']
results_dict = {}
for tissue in tissues:
    results = pd.read_csv(os.path.join("output/paired/rna_atac", tissue, "results.csv"), index_col=0)
    results["tissue"] = tissue
    results.index = tissue + "_" + results.index
    results_dict[tissue] = results
results = pd.concat(results_dict.values())
results

### Plotting the average default setting Alignment across datasets

### Averaging across all subsampled versions as well
default_results = results[results["var_num"]==2000]
default_results = default_results.drop(columns=["filter_method", "var_num", "cor_num", "cell_num", "index"])
default_results = default_results.groupby(["algorithm", "tissue"], as_index=False).mean()
default_results["tissue"] = pd.Categorical(default_results["tissue"], categories = tissues, ordered=True)
default_results

plt.figure(figsize=(12, 6))
sns.barplot(data=default_results, x='tissue', y='ari', hue='algorithm')

plt.title('Algorithm Performance by Tissue - Label Transfer(Higher is Better)')
plt.xlabel('Tissue')
plt.ylabel('Alignment Score')

plt.xticks(rotation=45, ha='right')

# Adjust the legend to be outside the plot
plt.legend(title='Algorithm', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()

### Plotting percent improvement of FOSCTTM  across tissues and conditions
# Filter and process default results
default_results = results[results["var_num"] == 2000].drop(columns=["filter_method", "index"])
default_results = default_results.groupby(["algorithm", "tissue", "var_num", "cor_num", "cell_num"], as_index=False).mean()
default_results["tissue"] = pd.Categorical(default_results["tissue"], categories=tissues, ordered=True)
# Filter and process filtered results
filtered_results = results[results["var_num"] != 2000].copy()
filtered_results["tissue"] = pd.Categorical(filtered_results["tissue"], categories=tissues, ordered=True)
filtered_results
# # Merge the two DataFrames
merged = filtered_results.merge(
    default_results[['algorithm', 'tissue', 'cell_num', 'foscttm', 'ari', 'asw', 'accuracy', 'asw-batch']],
    on=['algorithm', 'tissue', 'cell_num'],
    suffixes=('', '_default')
)
# # Calculate the percentage difference
# # Add the calculated column to filtered_results
# # Calculate raw improvements
metrics = ['foscttm', 'ari', 'asw', 'accuracy', 'asw-batch']
for metric in metrics:
    merged[metric+"_percent"] = (merged[metric] - merged[metric+'_default'])/merged[metric+'_default'].replace(0, np.nan)
merged


plt.figure(figsize=(20, 10))
sns.boxplot(data=merged, x='tissue', y='asw_percent', hue="var_num")
plt.title('Relative Improvement - ASW')

# Add a red dotted line at y=0
plt.axhline(y=0, color='red', linestyle='--', linewidth=1.5)
plt.axhline(y=0.1, color='red', linestyle='--', linewidth=1.5)
# plt.axhline(y=0.2, color='red', linestyle='--', linewidth=1.5)

plt.tight_layout()
plt.savefig("plots/benchmark/ASW_parameter_robustness.svg", format="svg", dpi=300)
plt.savefig("plots/benchmark/ASW_parameter_robustness.jpeg", format="jpeg", dpi=300)

