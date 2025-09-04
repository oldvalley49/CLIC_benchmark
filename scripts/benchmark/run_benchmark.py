import numpy as np
import pandas as pd
import os
import sys
import re
import scipy
import sklearn.metrics
from utility import *
# set current directory

# benchmark for paired data
output_dir = "output/paired/rna_atac"

args = sys.argv[1:]
tissue = args[0]
    
coembed_dir = os.path.join(output_dir,tissue, "coembed")
print(coembed_dir)
methods = [method for method in os.listdir(coembed_dir) if os.path.isfile(os.path.join(coembed_dir, method)) and not method.startswith(".")]
annotation_dir = os.path.join(output_dir, tissue, "ct_annotation")
if os.path.exists(os.path.join(output_dir, tissue, "results.csv")):
    result=pd.read_csv(os.path.join(output_dir, tissue, "results.csv"), index_col=0)
else:
    result = pd.DataFrame(columns=["algorithm", "filter_method", "var_num", "cor_num", "cell_num", "index", "foscttm", "ari", "asw", "accuracy", "asw-batch", "knn_auc"])
for method in methods:
    print(method)
    if method in result.index:
        continue
    if os.path.exists(os.path.join(output_dir, tissue, "ct_annotation", method))==False:
        continue
    # parse string
    temp = {}
    parse_result = re.split(r'[_\.]', method)
    temp["algorithm"] = parse_result[0]
    temp["filter_method"] = parse_result[1]
    temp["var_num"] = parse_result[2]
    temp["cor_num"] = parse_result[3]
    temp["cell_num"] = parse_result[4]
    temp["index"] = parse_result[5]
    rna, atac = load_coembed(os.path.join(output_dir,tissue, "coembed", method))
    coordinates_combined = pd.concat([rna, atac])
    annotation = pd.read_csv(os.path.join(output_dir, tissue, "ct_annotation", method), index_col=0)
    annotation_combined = pd.concat([annotation, annotation])
    rna_num = rna.shape[0]
    atac_num = atac.shape[0]
    batch_label = pd.Series(["RNA"] * rna_num + ["ATAC"] * atac_num)
    batch_label = pd.Series(["RNA"] * rna_num + ["ATAC"] * atac_num)
    
    index_combined = pd.Index(list(rna.index)+list(atac.index))
    annotation_combined.index = index_combined
    batch_label.index = index_combined

    temp["foscttm"] = foscttm(rna, atac)
    temp["ari"] = ari(annotation["ct_truth"], annotation["predicted.id"])
    temp["asw"] = asw(coordinates_combined.to_numpy(), annotation_combined["ct_truth"].values)
    temp["accuracy"] = accuracy(annotation["ct_truth"], annotation["predicted.id"])
    temp["asw-batch"] = asw_batch(coordinates_combined.to_numpy(), batch_label,annotation_combined["ct_truth"].values)
    temp["knn_auc"] = knn_auc(rna, atac)
    temp_df = pd.DataFrame(temp, index=[method])
    result = pd.concat([result, temp_df])

result.to_csv(os.path.join(output_dir, tissue, "results.csv"))