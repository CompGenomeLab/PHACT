import os
import sys


def make_query_fasta(acc_no,all_eu_file,output_path):

    protein_dict = {}

    with open(all_eu_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.split(">")[1].strip()
                if header.split("|")[1].split("|")[0] == acc_no:
                    line = next(f).strip()
                    protein_dict[header] = ''
                    while not line.startswith('>'):
                        protein_dict[header] += line
                        line = next(f).strip()
    f.close()

    with open(output_path,'w') as f:
        for header,sequence in protein_dict.items():
            f.write(">"+header+"\n" + sequence + "\n")
    f.close()


if __name__ == "__main__":
    acc_no = sys.argv[1]
    all_eu_file = sys.argv[2]
    output_path = sys.argv[3]
    make_query_fasta(acc_no,all_eu_file,output_path)
    
    



