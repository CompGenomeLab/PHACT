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
                                #print(best_hit)
                                blast_list.append(best_hit)
    return blast_list

def write_to_file(blast_list,all_eu_file,blastp_out_file,query_fasta):
    query_fasta_dict = get_fasta_dict(query_fasta)
    all_eu_dict = get_fasta_dict(all_eu_file)
    with open(blastp_out_file.split(".")[0]+".fasta",'w') as f:
        for k,v in all_eu_dict.items():
            #print(k.split(" ")[0])
            if k.split(" ")[0] in blast_list:
                f.write(">"+k + "\n" + v +"\n")
        for h,s in query_fasta_dict.items():
            f.write(">"+h + "\n" + s +"\n")
    f.close()


if __name__ == "__main__":
    blastp_out_file = sys.argv[1]
    blast_hit_number = sys.argv[2]
    max_e_value = sys.argv[3]
    min_identity = sys.argv[4]
    max_identity = sys.argv[5]
    blastdb_file = sys.argv[6]
    query_fasta = sys.argv[7]
    blast_list =  parse_blastp_out(blastp_out_file,blast_hit_number,max_e_value,min_identity,max_identity)
    write_to_file(blast_list,blastdb_file,blastp_out_file,query_fasta)
