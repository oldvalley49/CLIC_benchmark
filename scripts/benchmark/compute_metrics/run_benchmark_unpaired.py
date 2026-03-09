import numpy as np
import pandas as pd
import os
import sys
import re
import scipy
import sklearn.metrics
from utility import *


# benchmark for unpaired data
output_dir = "output/unpaired/rna_atac"

    
pairs = [pair for pair in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, pair)) and not pair.startswith(".")]
annotation_dir = os.path.join(output_dir, "ct_annotation")
for pair in pairs:

    if os.path.exists(os.path.join(output_dir, pair, "results.csv")):
        result=pd.read_csv(os.path.join(output_dir, pair, "results.csv"), index_col=0)
    else:
        result = pd.DataFrame(columns=["algorithm", "rna", "atac", "filter_method", "var_num", "cor_num", "rna_num", "atac_num", "batch", "ari", "asw", "asw-batch"])
    print(pair)
    methods = [i for i in os.listdir(os.path.join(output_dir, pair, "coembed")) if not i.startswith(".")]
    for method in methods:
        if method in result.index:
            continue
        if os.path.exists(os.path.join(output_dir, pair, "ct_annotation", method))==False:
            continue
        # parse string
        temp = {}
        parse_result = re.split(r'[_\.]', method)               
        temp["algorithm"] = parse_result[0]
        temp["rna"] = parse_result[1]
        temp["atac"] = parse_result[2]
        temp["filter_method"] = parse_result[3]
        temp["var_num"] = int(parse_result[4])
        temp["cor_num"] = parse_result[5]
        temp["rna_num"] = int(parse_result[6])
        temp["atac_num"] = int(parse_result[7])
        temp["batch"] = parse_result[8]
        rna, atac = load_coembed_unpaired(os.path.join(output_dir,pair,"coembed", method),temp["rna_num"])
        coordinates_combined = pd.concat([rna, atac])
        annotation_atac = pd.read_csv(os.path.join(output_dir,pair, "ct_annotation", method), index_col=0)

        # get rna cells annotation
        annotation_rna = pd.read_csv(os.path.join("data/processed_data/rna_atac/", temp["rna"], "ct_annotation.csv"), index_col=0)
        barcodes_rna = pd.read_csv(os.path.join("data/processed_data/rna_atac/", temp["rna"], "barcode.csv"), index_col=0)
        annotation_rna.index = barcodes_rna["x"]
        rna.index = rna.index.str.replace(r'^(rna_|RNA_)', '', regex=True)
        assert set(rna.index) == set(annotation_rna.index), "Mismatch between indices!"
        annotation_rna = annotation_rna.loc[rna.index]

        annotation_combined = pd.concat([annotation_rna["x"], annotation_atac["ct_truth"]])

        rna_num = rna.shape[0]
        atac_num = atac.shape[0]
        batch_label = pd.Series(["RNA"] * rna_num + ["ATAC"] * atac_num)
        batch_label = pd.Series(["RNA"] * rna_num + ["ATAC"] * atac_num)
        batch_label.index = coordinates_combined.index

        temp["ari"] = ari(annotation_atac["ct_truth"], annotation_atac["predicted.id"])
        temp["asw"] = asw(coordinates_combined.to_numpy(), annotation_combined.values)
        temp["asw-batch"] = asw_batch(coordinates_combined.to_numpy(), batch_label, annotation_combined.values)
        temp_df = pd.DataFrame(temp, index=[method])
        result = pd.concat([result, temp_df])           

    result.to_csv(os.path.join(output_dir, pair, "results.csv"))