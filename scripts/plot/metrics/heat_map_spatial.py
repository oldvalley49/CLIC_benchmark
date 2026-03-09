import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
import numpy as np

# load benchmarking results for each run
tissues = ["brain", "Melanoma", "mEmbryo", "mEmbryo2", "mBrain2"]
results_dict = {}
for tissue in tissues:
    file_path = os.path.join("output/spatial", tissue, "results.csv")
    if os.path.exists(file_path):
        results = pd.read_csv(file_path, index_col=0)
    else:
        continue
    results["tissue"] = tissue
    results.index = tissue + "_" + results.index
    results_dict[tissue] = results
results = pd.concat(results_dict.values())
results


### Plotting improvement across datasets

# Filter and process default results
default_results = results[(results["var_num"] == 2000) & (results["filter_method"] == "corr")].drop(columns=["filter_method", "index"])
default_results = default_results.groupby(["algorithm", "tissue", "var_num", "cor_num", "cell_num"], as_index=False).mean()
default_results["tissue"] = pd.Categorical(default_results["tissue"], categories=tissues, ordered=True)
# Filter and process filtered results
filtered_results = results[(results["var_num"] >= 4000) & (results["var_num"]<= 8000)].copy()
filtered_results["tissue"] = pd.Categorical(filtered_results["tissue"], categories=tissues, ordered=True)
filtered_results
# Merge the two DataFrames
merged = filtered_results.merge(
    default_results[['algorithm', 'tissue', 'cell_num', 'foscttm','knn_auc']],
    on=['algorithm', 'tissue', 'cell_num'],
    suffixes=('', '_default')
)
# Compute Improvement
merged['foscttm_improvement'] = (merged['foscttm']- merged['foscttm_default'])

metrics = ['foscttm','knn_auc']
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
    fig, ax = plt.subplots(figsize=(10, 10))

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
    plt.savefig(f"plots/benchmark/spatial_heatmap_corr_{metric}.svg")
    plt.savefig(f"plots/benchmark/spatial_heatmap_corr_{metric}.png")
    plt.close()






