import os
import sys
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

if __name__ == "__main__":
    fasta_file  = sys.argv[1]

