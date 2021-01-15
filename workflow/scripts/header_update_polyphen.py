import os
import sys
import re

def read_alignment(aln_file):
    seq_dict = {}
    i=0
    with open(aln_file) as f:
        for line in f:
            if line.startswith("sp|") or line.startswith("tr|"):
                header = line.split("  ")[0].strip()
                header = re.sub(r'[^\w>]', '_',header)
                seq = ''.join(line.split("  ")[1:])
                seq = seq.lstrip()
                seq = seq.strip()
                seq_dict[header] = seq
    return seq_dict

def write_new_aln(seq_dict,new_file):
    with open(new_file,'w') as f:
        for k,v in seq_dict.items():
            f.write(">" + k + "\n" + v + "\n")
    f.close()

if __name__ == "__main__":
    aln_file = sys.argv[1]
    seq_dict = read_alignment(aln_file)
    new_file = sys.argv[2]
    write_new_aln(seq_dict,new_file)