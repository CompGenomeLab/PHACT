import sys
import json
import os
import re
from fasta_dict import *

def change_name(fasta_file):
    fasta_dict1 = get_fasta_dict(fasta_file)
    new_dict = {}
    for a,v in fasta_dict1.items():
        new_a = a.split(" ")[0]+"_"+a.split("OX=")[1].split(" ")[0]
        new_a = re.sub(r'[^\w>]', '_',new_a)
        new_dict[new_a] = v

    with open(fasta_file.split(".")[0]+"_new_header.fasta",'w') as f:
        for key, value in new_dict.items():
            f.write(">"+key+"\n"+value+"\n")

    f.close()



if __name__ == "__main__":
    fasta_file = sys.argv[1]
    change_name(fasta_file)

