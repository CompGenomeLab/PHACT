import os
import sys
import logging
import subprocess
import csv
import json
import tempfile
import logging
import threading

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



def get_acc_numbers(uniprot_ids_file,base_path):
    acc_no_list = []
    with open(base_path + uniprot_ids_file,'r') as f:
        lines = f.readlines()
        for line in lines:
            acc_no_list.append(line.strip())
    return acc_no_list

def make_directories(acc_no_list):
    for acc_no in acc_no_list:
        if not os.path.isdir(acc_no):
            os.makedirs(acc_no)
    
    
def get_fasta_files(acc_no,protein_directory,all_eu_dict):
    protein_dict = {}
    for key, value in all_eu_dict.items():
        key_acc_no = key.split("|")[1].split("|")[0]
        if acc_no == key_acc_no:
            protein_dict[key] = value
    with open(protein_directory + acc_no+".fasta",'w') as f:
        for header,sequence in protein_dict.items():
            f.write(">"+header+"\n" + sequence + "\n")
    f.close()

def write_config_file(acc_no,base_path,weights):
    protein_directory = base_path + acc_no
    with open(protein_directory+"/config.yml",'w') as f:
        config= 'configfile: "config.yml"'
        channels = "channels: "
        blast = "  -blast"
        python = "  -python"
        dependencies = "dependencies: "
        blast_= "  - blast = 2.9.0"
        python_ = "  - python = 3.7.4"
        text = "workdir: "  + '"' +  protein_directory + '"'
        output_name = "output_name: " + '"'+  acc_no + '"'
        fasta_file = "fasta_file: " +'"'+ acc_no + ".fasta" + '"'
        fasta_folder = "fasta_folder: " + '"' + base_path + acc_no + '"'
        blastdb_folder = 'blastdb_folder: "/cta/users/abircan/blastdb"'
        blastdb_file = 'blastdb_file: "all_eu.fasta"'
        aa_rate_file_folder= 'aa_rate_file_folder: "/cta/users/abircan"'
        aa_rate_file = 'aa_rate_file: "JTT.dat.txt"'
        weight = 'weights: ' +'"'+ weights + '"'
        pattern = "pattern: "+ '"'+ weights[-1]+ '"'
        L = [config,"\n",channels,"\n",blast,"\n",python,"\n",dependencies,"\n",blast_,"\n",python_,"\n", text,"\n",output_name,"\n", fasta_file, "\n",fasta_folder,"\n",blastdb_folder,"\n",blastdb_file,"\n",aa_rate_file_folder,"\n",aa_rate_file,"\n",weight,"\n",pattern,"\n"]
        f.writelines(L) 
    f.close()


def create_snakemake_folder(acc_no_list,base_path,weights):
    for acc_no in acc_no_list:
        protein_directory = base_path + acc_no +"/"
        get_fasta_files(acc_no,protein_directory,all_eu_dict)
        write_config_file(acc_no,base_path,weights)



if __name__ == "__main__":
    uniprot_ids_file = sys.argv[1]
    base_path = sys.argv[2]
    weights = sys.argv[3]
    all_eu_dict = get_fasta_dict("/cta/users/abircan/blastdb/all_eu.fasta")
    acc_no_list = get_acc_numbers(uniprot_ids_file,base_path)
    make_directories(acc_no_list)
    create_snakemake_folder(acc_no_list,base_path,weights)


