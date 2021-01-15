rule run_codeml:
    input:
        no_gap_file = "{output_folder}/results/{msa_name}/2_msa/{msa_name}_nogap_msa.fasta",
        tree = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.bestTree_unrooted",
    params:
        wdir = "{output_folder}/results/{msa_name}/4_codeml",
    output:
        codeml_out="{output_folder}/results/{msa_name}/4_codeml/{msa_name}_codeml",
    conda:
        "../envs/paml.yml"
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_run_codeml.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_run_codeml.out"
    resources:
        time=1200
    shell:
        "PATH={config[output_folder]}/resources/paml4.9j/bin:$PATH && python scripts/codeml.py {input.no_gap_file} {input.tree} {output.codeml_out} {params.wdir} {config[seqtype]} {config[verbose]} {config[noisy]} {config[clock]} {config[aa_dist]} {config[aa_rate_file]} {config[model]} {config[icode]} {config[Mgene]} {config[fix_alpha]} {config[alpha]} {config[Malpha]} {config[ncatG]} {config[getSE]} {config[RateAncestor]} {config[Small_Diff]} {config[cleandata]} {config[method]} 2> {log}"
