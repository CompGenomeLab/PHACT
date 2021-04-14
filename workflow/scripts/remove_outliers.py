from ete3 import Tree
import os
import sys
from fasta_dict import *
import scipy
import pandas as pd
import numpy as np
from scipy.stats import zscore
import scipy.stats as stat
from scipy.stats import iqr

def get_leaf_lenghts(tree_file):
    t = Tree(tree_file,format=1)
    leaf_lengths = [leaf.dist for leaf in t.get_leaves()]
    print(type(leaf_lengths))
    return leaf_lengths

def get_outliers(leaf_lengths,max_deviations):
    df = np.array(leaf_lengths)
    mean = np.mean(df)
    standard_deviation = np.std(df)
    distance_from_mean = abs(df - mean)
    not_outlier = distance_from_mean < max_deviations * standard_deviation
    no_outliers = df[not_outlier]
    no_outliers_list = no_outliers.tolist()
    #print(no_outliers_list)
    return no_outliers_list

def remove_outliers(tree_file,msa_file,no_outliers_list,no_outlier_tree,no_outlier_msa):
    t = Tree(tree_file,format=1)
    msa_dict = get_fasta_dict(msa_file)
    prune_list = []
    for leaf in t:
        if leaf.dist in no_outliers_list:
            prune_list.append(leaf.name)

    t.prune(prune_list)
    t.write(outfile = no_outlier_tree,format=1)

    with open(no_outlier_msa,'w') as f:
        for header, seq in msa_dict.items():
            if header in prune_list:
                f.write(">" + header + "\n" + seq + "\n")
    f.close()
    

if __name__ == "__main__":
    tree_file = sys.argv[1]
    msa_file = sys.argv[2]
    no_outlier_tree = sys.argv[3]
    no_outlier_msa = sys.argv[4]
    max_deviations = int(sys.argv[5])
    leaf_lengths = get_leaf_lenghts(tree_file)
    no_outliers_list = get_outliers(leaf_lengths,max_deviations)
    remove_outliers(tree_file,msa_file,no_outliers_list,no_outlier_tree,no_outlier_msa)





