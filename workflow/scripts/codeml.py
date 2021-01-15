from Bio.Phylo.PAML import codeml
import sys
import os

def write_tcl_file(alignment_file,tree_file,out_file,directory,seqtype,verbose,noisy,clock,aaDist,aa_rate_file,model,icode,Mgene,fix_alpha, alpha,Malpha,ncatG,getSE,RateAncestor,Small_Diff,cleandata,method):
    tcl = codeml.Codeml()
    tcl.alignment = alignment_file
    tcl.tree = tree_file
    tcl.out_file = out_file
    tcl.working_dir = directory
    tcl.set_options(seqtype=int(seqtype),
            verbose=int(verbose),
	        noisy=int(noisy),
            clock = int(clock),
            aaDist = int(aaDist),
            aaRatefile = aa_rate_file,
            model = int(model),
            icode = int(icode),
            Mgene = int(Mgene),
            fix_alpha = int(fix_alpha),   
            alpha = float(alpha), 
            Malpha = int(Malpha),   
            ncatG = int(ncatG),  
            getSE = int(getSE),  
            RateAncestor = int(RateAncestor),  
            Small_Diff = float(Small_Diff),
            cleandata = int(cleandata),  
            method = int(method)) 
    tcl.run() 

if __name__ == "__main__":
    alignment_file = sys.argv[1]
    tree_file = sys.argv[2]
    out_file = sys.argv[3]
    directory = sys.argv[4]
    seqtype = sys.argv[5]
    verbose = sys.argv[6]
    noisy = sys.argv[7]
    clock = sys.argv[8]
    aaDist = sys.argv[9]
    aa_rate_file = sys.argv[10]
    model = sys.argv[11]
    icode = sys.argv[12]
    Mgene = sys.argv[13]
    fix_alpha = sys.argv[14]
    alpha = sys.argv[15] 
    Malpha = sys.argv[16]   
    ncatG = sys.argv[17]  
    getSE = sys.argv[18]  
    RateAncestor = sys.argv[19]
    Small_Diff = sys.argv[20]
    cleandata = sys.argv[21] 
    method = sys.argv[22]
    write_tcl_file(alignment_file,tree_file,out_file,directory,seqtype,verbose,noisy,clock,aaDist,aa_rate_file,model,icode,Mgene,fix_alpha, alpha,Malpha,ncatG,getSE,RateAncestor,Small_Diff,cleandata,method)





