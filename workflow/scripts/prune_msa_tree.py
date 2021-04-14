import os
import sys
import ete3
from ete3 import Tree
import re

def get_query_id(fasta_file):
    with open(str(fasta_file),'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.split(">")[1].strip()
                query_id = header.split(" ")[0]+"_"+header.split("OX=")[1].split(" ")[0]
                query_id = re.sub(r'[^\w>]', '_',query_id)
    #print(query_id)
    return query_id


def get_fasta_dict(fasta_file):
    fasta_dict = {}
    
    with open(fasta_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.split(">")[1].strip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.strip() 
    f.close()
    return fasta_dict

def parse_blastp_out(blastp_out_file,blast_hit_number,max_e_value,min_identity,max_identity):
    blast_list = []
    with open(blastp_out_file,'r') as filein_:
        for line in filein_:
            if len(blast_list) <= int(blast_hit_number):
                if line.startswith(">"):
                    best_hit = line.split(">")[1].strip()
                    best_hit = best_hit.split(" ")[0]
                    #print(best_hit)
                    while not line.startswith(" Score"):
                        line = next(filein_)
                        if line.startswith(" Score"):
                            #print("???")
                            score = float(line.split("Score = ")[1].split(" bits")[0])
                            e_value = float(line.split("Expect = ")[1].split(",")[0])
                            identity_line = next(filein_)
                            identity = float(identity_line.split("Identities = ")[1].split("(")[1].split("%")[0])
                            if e_value < float(max_e_value) and float(max_identity)>identity>float(min_identity):
                                best_hit = best_hit.split("|")[1].split("|")[0]
                                blast_list.append(best_hit)
    #print(blast_list)
    return blast_list

def prune_msa_tree(ml_tree,pruned_tree,msa_file,pruned_msa,blast_list,query_id):
    t = Tree(ml_tree,format=1)
    msa_dict = get_fasta_dict(msa_file)
    pruned_dict = {}
    prune_list = []
    for leaf in t:
        if leaf.name.split("_")[1].split("_")[0] in blast_list:
            prune_list.append(leaf.name)
    if not query_id in prune_list:
        prune_list.append(query_id)
    t.prune(prune_list)
    t.write(outfile=pruned_tree,format=1)

    for k,v in msa_dict.items():
        if k.split("_")[1].split("_")[0] in blast_list:
            pruned_dict[k] = v
    if not query_id in pruned_dict.keys():
        pruned_dict[query_id] = msa_dict[query_id]
    
    with open(pruned_msa,'w') as f:
        for k,v in pruned_dict.items():
            f.write(">" + k + "\n" + v + "\n")
    f.close()

if __name__ == "__main__":
    blastp_out_file = sys.argv[1]
    blast_hit_number = sys.argv[2]
    max_e_value = sys.argv[3]
    min_identity = sys.argv[4]
    max_identity = sys.argv[5]
    ml_tree = sys.argv[6]
    pruned_tree = sys.argv[7]
    msa_file = sys.argv[8]
    pruned_msa = sys.argv[9]
    fasta_file = sys.argv[10]
    query_id = get_query_id(fasta_file)
    blast_list = parse_blastp_out(blastp_out_file,blast_hit_number,max_e_value,min_identity,max_identity)
    prune_msa_tree(ml_tree,pruned_tree,msa_file,pruned_msa,blast_list,query_id)





        
    





