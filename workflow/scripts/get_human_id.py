import re
import sys

def get_human_id(fasta_file):
    with open(str(fasta_file),'r') as f:
        for line in f:
            if line.startswith('>'):
                line = line[0] + line[1:].replace('>', '-') 
                header = line.split(">")[1].strip()
                human_id = header.split(" ")[0]+"_"+header.split("OX=")[1].split(" ")[0]
                human_id = re.sub(r'[^\w>]', '_',human_id)
    print(human_id)
    return human_id

get_human_id(sys.argv[1])