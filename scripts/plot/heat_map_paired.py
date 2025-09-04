import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
import numpy as np

# load benchmarking results for each run
tissues = ["PBMC", "BMMC-s1d2", "BMMC-s4d8", "TDBM", "mBrain", 'mSkin', 'mRetina']
results_dict = {}
for tissue in tissues:
    file_path = os.path.join("output/paired/rna_atac", tissue, "results.csv")
    if os.path.exists(file_path):
        results = pd.read_csv(file_path, index_col=0)
    else:
        continue
    results["tissue"] = tissue
    results.index = tissue + "_" + results.index
    results_dict[tissue] = results
results = pd.concat(results_dict.values())
results

### Plotting the average default setting Alignment across datasets

### Averaging across all subsampled versions as well
default_results = results[(results["var_num"] == 2000) & (results["filter_method"] == "corr")]
default_results = default_results.drop(columns=["filter_method", "var_num", "cor_num", "cell_num", "index"])
default_results = default_results.groupby(["algorithm", "tissue"], as_index=False).mean()
default_results["tissue"] = pd.Categorical(default_results["tissue"], categories = tissues, ordered=True)
default_results

plt.figure(figsize=(12, 6))
sns.barplot(data=default_results, x='tissue', y='foscttm', hue='algorithm')

plt.title('Algorithm Performance by Tissue - FOSCTTM (Higher is Better)')
plt.xlabel('Tissue')
plt.ylabel('FOSCTTM')

plt.xticks(rotation=45, ha='right')
plt.legend(title='Algorithm', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("plots/benchmark/avg_default_sc.png")


### Plotting improvement across datasets

# Filter and process default results
default_results = results[(results["var_num"] == 2000) & (results["filter_method"] == "corr")].drop(columns=["filter_method", "index"])
default_results = default_results.groupby(["algorithm", "tissue", "var_num", "cor_num", "cell_num"], as_index=False).mean()
default_results["tissue"] = pd.Categorical(default_results["tissue"], categories=tissues, ordered=True)
# Filter and process filtered results
filtered_results = results[results["var_num"] >= 4000].copy()
filtered_results["tissue"] = pd.Categorical(filtered_results["tissue"], categories=tissues, ordered=True)
filtered_results
# Merge the two DataFrames
merged = filtered_results.merge(
    default_results[['algorithm', 'tissue', 'cell_num', 'foscttm', 'ari', 'asw', 'accuracy', 'asw-batch', 'knn_auc']],
    on=['algorithm', 'tissue', 'cell_num'],
    suffixes=('', '_default')
)
# Compute Improvement
merged['foscttm_improvement'] = (merged['foscttm']- merged['foscttm_default'])

metrics = ['foscttm', 'ari', 'asw', 'accuracy', 'asw-batch', 'knn_auc']
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
    plt.savefig(f"plots/benchmark/paired_heatmap_corr_{metric}.svg")
    plt.close()




# print corr and single side by side
for metric in metrics:
    plot_metric = metric+"_improvement"
    # Filter the data for 'single' and 'corr' filter_method
    corr_data = merged[merged['filter_method'] == 'single']
    single_data = merged[merged['filter_method'] == 'corr']

    # Group and pivot for corr
    grouped_corr = corr_data.groupby(['varnum_method', 'tissue'])[plot_metric].mean().reset_index()
    pivot_corr = grouped_corr.pivot(index='varnum_method', columns='tissue', values=plot_metric)
    pivot_corr = pivot_corr.reindex(sorted(pivot_corr.index, key=lambda x: (int(x.split("_")[0]), x.split("_")[1])))
    pivot_corr = pivot_corr.reindex(columns=tissues)

    # Group and pivot for single
    grouped_single = single_data.groupby(['varnum_method', 'tissue'])[plot_metric].mean().reset_index()
    pivot_single = grouped_single.pivot(index='varnum_method', columns='tissue', values=plot_metric)
    pivot_single = pivot_single.reindex(sorted(pivot_single.index, key=lambda x: (int(x.split("_")[0]), x.split("_")[1])))
    pivot_single = pivot_single.reindex(columns=tissues)

    # Define the normalization
    vmin = min(pivot_corr.min().min(), pivot_single.min().min())
    vmax = max(pivot_corr.max().max(), pivot_single.max().max())
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    # Create a figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    # Heatmap for corr
    sns.heatmap(
        pivot_corr,
        annot=True,
        fmt=".4f",
        cmap="coolwarm",  # Diverging colormap
        norm=norm,        # Normalize to center 0
        cbar_kws={'label': plot_metric},
        ax=axes[0]
    )
    axes[0].set_title(plot_metric + " (corr)")
    axes[0].set_xlabel("Tissue")
    axes[0].set_ylabel("Var_Num and Algorithm")
    axes[0].tick_params(axis='x', rotation=45)

    # Heatmap for single
    sns.heatmap(
        pivot_single,
        annot=True,
        fmt=".4f",
        cmap="coolwarm",  # Diverging colormap
        norm=norm,        # Normalize to center 0
        cbar_kws={'label': plot_metric},
        ax=axes[1]
    )
    axes[1].set_title(plot_metric + " (single)")
    axes[1].set_xlabel("Tissue")
    axes[1].set_ylabel("Var_Num and Algorithm")
    axes[1].tick_params(axis='x', rotation=45)

    # Adjust layout
    plt.tight_layout()
    plt.savefig(f"plots/benchmark/paired_heatmap_{metric}.png")

