import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
import numpy as np

# load benchmarking results for each run
tissues_rna = ["BMMC-s1d2", "BMMC-s4d8", "BMMC-s3d6"]
tissues_atac = ["BMMC-s2d4", "BMMC-s1d2", "BMMC-s2d4"]

tissues_zipped = list(zip(tissues_rna, tissues_atac))
tissues = list()
results_dict = {}
for rna, atac in tissues_zipped:
    pair = rna + "+" + atac
    tissues.append(pair)
    file_path = os.path.join("output/unpaired/rna_atac", pair, "results.csv")
    if os.path.exists(file_path):
        results = pd.read_csv(file_path, index_col=0)
    else:
        continue
    results["tissue"] = pair
    results.index = pair + "_" + results.index
    results_dict[pair] = results
results = pd.concat(results_dict.values())
results

### Plotting the average default setting Alignment across datasets

### Averaging across all subsampled versions as well
# default_results = results[(results["var_num"] == 2000) & (results["filter_method"] == "corr")]
# default_results = default_results.drop(columns=["filter_method", "var_num", "cor_num", "rna_num", "atac_num", "batch", "rna", "atac"])
# default_results = default_results.groupby(["algorithm", "tissue"], as_index=False).mean()
# default_results["tissue"] = pd.Categorical(default_results["tissue"], categories = tissues, ordered=True)
# default_results

# plt.figure(figsize=(12, 6))
# sns.barplot(data=default_results, x='tissue', y='foscttm', hue='algorithm')

# plt.title('Algorithm Performance by Tissue - FOSCTTM (Higher is Better)')
# plt.xlabel('Tissue')
# plt.ylabel('FOSCTTM')

# plt.xticks(rotation=45, ha='right')
# plt.legend(title='Algorithm', bbox_to_anchor=(1.05, 1), loc='upper left')

# plt.tight_layout()
# plt.savefig("plots/benchmark/avg_default_sc.svg")


### Plotting improvement across datasets

# Filter and process default results
default_results = results[(results["var_num"] == 2000) & (results["filter_method"] == "corr")].drop(columns=["filter_method", "batch"])
default_results = default_results.groupby(["algorithm", "tissue", "var_num", "cor_num", "rna_num", "atac_num", "rna", "atac"], as_index=False).mean()
default_results["tissue"] = pd.Categorical(default_results["tissue"], categories=tissues, ordered=True)
# Filter and process filtered results
filtered_results = results[results["var_num"] >= 4000].copy()
filtered_results["tissue"] = pd.Categorical(filtered_results["tissue"], categories=tissues, ordered=True)
filtered_results
# Merge the two DataFrames
merged = filtered_results.merge(
    default_results[['algorithm', 'tissue', 'rna_num', 'atac_num', 'ari', 'asw', 'asw-batch']],
    on=['algorithm', 'tissue', 'rna_num', 'atac_num'],
    suffixes=('', '_default')
)

metrics = ['ari', 'asw', 'asw-batch']
for metric in metrics:
    merged[metric+"_improvement"] = (merged[metric] - merged[metric+'_default'])
merged['varnum_method'] = merged['var_num'].astype(str) + "_" + merged['algorithm']

for metric in metrics:
    plot_metric = metric + "_improvement"   
    corr_data = merged[merged['filter_method'] == 'corr']

    # Group and pivot for corr
    grouped_corr = corr_data.groupby(['varnum_method', 'tissue'])[plot_metric].mean().reset_index()
    pivot_corr = grouped_corr.pivot(index='varnum_method', columns='tissue', values=plot_metric)
    pivot_corr = pivot_corr.reindex(sorted(pivot_corr.index, key=lambda x: (int(x.split("_")[0]), x.split("_")[1])))
    pivot_corr = pivot_corr.reindex(columns=tissues)

    # Define normalization
    vmin = pivot_corr.min().min()
    vmax = pivot_corr.max().max()
    max_num = max(abs(vmin), abs(vmax))
    norm = TwoSlopeNorm(vmin=-max_num, vcenter=0, vmax=max_num)

    # Create single subplot
    fig, ax = plt.subplots(figsize=(6, 10))

    # Heatmap for corr
    sns.heatmap(
    pivot_corr,
    annot=True,
    fmt=".4f",
    cmap="coolwarm",
    norm=norm,
    cbar_kws={'label': plot_metric},
    ax=ax
    )
    ax.set_title(f"{plot_metric} (corr)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Parameter-Algorithm Combination")
    ax.tick_params(axis='x', rotation=45)

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(f"plots/benchmark/BMMCBatch_heatmap_corr_{metric}.svg")




