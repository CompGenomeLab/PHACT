from ete3 import Tree
import os
import sys
from fasta_dict import *
import re
import ete3
from ete3 import Tree

def unroot_tree(tree_file,no_gap_msa_file):
    msa_dict = get_fasta_dict(no_gap_msa_file)
    seqs_list = list(msa_dict.keys())
    t = Tree(tree_file,format=1)
    prune_list = []
    leaves = t.get_leaf_names()
    for leaf in leaves:
        if leaf in seqs_list:
            prune_list.append(leaf)
    
    t.prune(prune_list)
    t.unroot()
    t.write(outfile=tree_file+"_unrooted",format=1)

if __name__ == "__main__":
    tree_file = sys.argv[1]
    no_gap_msa_file = sys.argv[2]
    unroot_tree(tree_file,no_gap_msa_file)

