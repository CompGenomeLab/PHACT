import sys
import ete3 
from ete3 import Tree, PhyloTree
from fasta_dict import *
from scipy import stats
import numpy as np
from scipy.stats import ttest_ind, ttest_ind_from_stats


def remove_paralogs(tree_file,fasta_file,query_name,output_tree,output_fasta):
    t = Tree(tree_file,format=1)
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    for node in t.traverse(strategy="postorder"):
        if not node.is_leaf():
            children = node.get_children()
            if len(children) == 2:
                ch1,ch2 = children
                if ch1.is_leaf() and ch2.is_leaf():
                    if ch1.name.split("_")[-1] == ch2.name.split("_")[-1]:
                        if ch1.name.split("_")[-1] == "9606" or ch2.name.split("_")[-1] == "9606":
                            node.name = "HUMAN node"
                        if ch1.name == query_name or ch2.name == query_name:
                            node.name = "Query node"
                        else:
                            node.name = "Duplication node"

                else:
                    grand_leaves1_names = [leaf.name for leaf in ch1.get_leaves()]
                    grand_leaves2_names = [leaf.name for leaf in ch2.get_leaves()]
                    grand_leaves1 = set([leaf.name.split("_")[-1] for leaf in ch1.get_leaves()])
                    grand_leaves2 = set([leaf.name.split("_")[-1] for leaf in ch2.get_leaves()])
                    if query_name in grand_leaves1_names or query_name in grand_leaves2_names:
                            node.name = "Query node"
                    elif "9606" in grand_leaves1 or "9606" in grand_leaves2:
                        node.name = "HUMAN node"
                    elif grand_leaves1.intersection(grand_leaves2):
                        node.name = "Duplication node"

                
    new_t = t
    # #print(new_t)
    for node in new_t.traverse(strategy="postorder"):
        if node.name == "HUMAN node":
            node.detach()
        elif node.name == "Duplication node":
            children = node.get_children()
            if len(children) == 2:
                ch1, ch2 = children
                if ch1.is_leaf() and ch2.is_leaf():
                    ancestor = new_t.get_common_ancestor(ch1,ch2)
                    farthest = list(ancestor.get_farthest_leaf())[0]
                    node.remove_child(farthest)
                else:
                    grand_leaves1 = [leaf for leaf in ch1.get_leaves()]
                    grand_leaves2 = [leaf for leaf in ch2.get_leaves()]
                    grand_leaves1_distances = []
                    grand_leaves2_distances = []
                    for grand_leaf in grand_leaves1:
                        grand_leaves1_distances.append(node.get_distance(grand_leaf))    
                    for grand_leaf2 in grand_leaves2:
                        grand_leaves2_distances.append(node.get_distance(grand_leaf2))   
                    t, p = ttest_ind(grand_leaves1_distances, grand_leaves2_distances, equal_var=False)
                    avg1 = np.mean(grand_leaves1_distances)
                    avg2 = np.mean(grand_leaves2_distances)
                    if p < 0.0001:
                        if avg1 < avg2:
                            ch2.detach()
                        else:
                            ch1.detach()
                

    final_t= new_t
    old_fasta_dict = get_fasta_dict(fasta_file)
    new_fasta_dict = {}
    for node in final_t.traverse():
        leaf_names= node.get_leaf_names()
        for leaf in leaf_names:
            if leaf in old_fasta_dict.keys():
                new_fasta_dict[leaf] = old_fasta_dict[leaf]

    with open(output_fasta,'w') as f:
        for k,v in new_fasta_dict.items():
            f.write(">" + k + "\n" + v + "\n")
    f.close()

    prune_t = Tree(tree_file,format=1)
    final_fasta_dict = get_fasta_dict(output_fasta)
    prune_list = []
    for final_node in prune_t.traverse():
        final_leaf_names = final_node.get_leaf_names()
        for final_leaf in final_leaf_names:
            if final_leaf in final_fasta_dict.keys():
                if not final_leaf in prune_list:
                    prune_list.append(final_leaf)
    
    prune_t.prune(prune_list)
    new_node = prune_t.search_nodes()[0]
    new_leaf_names = new_node.get_leaf_names()
    new_node.write(outfile=output_tree,format=1)



if __name__ == "__main__":
    tree_file = sys.argv[1]
    fasta_file = sys.argv[2]
    query_name = sys.argv[3]
    output_tree = sys.argv[4]
    output_fasta = sys.argv[5]
    remove_paralogs(tree_file,fasta_file,query_name,output_tree,output_fasta)
