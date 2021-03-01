import os
import sys
from fasta_dict import *
import re

def get_human_id(query_file):
    query_dict = get_fasta_dict(query_file)
    for k,v in query_dict.items():
        human_id = k.split(" ")[0]+"_"+k.split("OX=")[1].split(" ")[0]
        human_id = re.sub(r'[^\w>]', '_',human_id)
    #print(human_id)
    return human_id

def get_gap_positions(seqDict,human_id):

    for key, value in seqDict.items():
        if human_id  in key:
            #print(value)
            indices = [i for i, x in enumerate(value) if x == "-"]
            #print(indices)
    return indices

def get_sequences_without_gap(seqDict,indices):
    new_seqDict = {}
    for k,v in seqDict.items():
        new_seqDict[k] = ''.join([v[i] for i in range(len(v)) if i not in indices])
        #print(new_seqDict)  
    return new_seqDict




def write_new_fasta(new_seqDict,fasta_file,output_file):
    #with open (fasta_file.split("_blasthits_new_header_msa")[0]+"_nogap_msa.fasta",'w') as new_file:
    with open (output_file,'w') as new_file:
        #print(new_file)
        for newheader,value in  new_seqDict.items():
            new_value = re.sub(r'[BXJZ]', '-',value)
            new_file.write(">" +newheader+"\n" + new_value+ "\n")
    new_file.close()



if __name__ == "__main__":
    query_file = sys.argv[1]
    fasta_file= sys.argv[2]
    output_file= sys.argv[3]
    human_id= get_human_id(query_file)
    seqDict = get_fasta_dict(fasta_file)
    gap_indices = get_gap_positions(seqDict,human_id)
    new_seqDict = get_sequences_without_gap(seqDict,gap_indices)
    write_new_fasta(new_seqDict,fasta_file,output_file)

