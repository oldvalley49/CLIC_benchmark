import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
import numpy as np

tissues = [i for i in os.listdir('output/paired/rna_atac') if not i.startswith('.')]
results_dict = {}
for tissue in tissues:
    results = pd.read_csv(os.path.join("output/paired/rna_atac", tissue, "results.csv"), index_col=0)
    results["tissue"] = tissue
    results.index = tissue + "_" + results.index
    results_dict[tissue] = results
results = pd.concat(results_dict.values())
results = results[results['filter_method'] == 'corr'].copy()

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

all_var_nums = sorted(merged['var_num'].unique())

algorithms = ['seurat', 'liger', 'bindsc', 'glue']
tissue_order = ['PBMC', 'BMMC-s1d2', 'BMMC-s4d8', 'TDBM', 'mBrain', 'mSkin', 'mRetina']

palette = sns.color_palette("ch:s=-.2,r=.6", n_colors=len(all_var_nums))


for algorithm in algorithms:
    merged_algorithm = merged[merged['algorithm'] == algorithm]
    # ensure var_num is categorical with all levels
    merged_algorithm['var_num'] = pd.Categorical(
        merged_algorithm['var_num'],
        categories=all_var_nums,
        ordered=True
    )
    plt.figure(figsize=(20, 10))
    sns.barplot(data=merged_algorithm, x='tissue', y='foscttm_percent', hue="var_num", order = tissue_order, palette = palette)
    plt.title(f'Relative Improvement - FOSCTTM - {algorithm}')
    plt.axhline(y=0, color='red', linestyle='--', linewidth=1.5)
    plt.tight_layout()
    plt.savefig(f"plots/benchmark/foscttm_parameter_robustness_{algorithm}.svg", format="svg", dpi=600)
    plt.savefig(f"plots/benchmark/foscttm_parameter_robustness_{algorithm}.jpeg", format="jpeg", dpi=600)

# find optimal parameter value

for algorithm in algorithms:
    merged_algorithm = merged[merged['algorithm'] == algorithm]
    df = merged_algorithm.copy()
    df_avg = (
        df.groupby(['tissue', 'var_num'], as_index=False)['foscttm_percent']
        .mean()
        .rename(columns={'foscttm_percent': 'foscttm_percent_mean'})
    )
    pt = df_avg.pivot(index='tissue', columns='var_num', values='foscttm_percent_mean')
    ranks = pt.rank(axis=1, ascending=False, method='average')
    k = pt.notna().sum(axis=1)
    for t in ranks.index:
        ranks.loc[t, ranks.loc[t].isna()] = k.loc[t]
    inv_scores = (k - ranks + 1).div(k, axis=0)
    final_scores = ranks.mean(axis=0).to_frame('avg_rank_score')
    final_ranking = final_scores.sort_values('avg_rank_score', ascending=True)
    print(f"parameter ranking, {algorithm}")
    print(final_ranking.head(10))
