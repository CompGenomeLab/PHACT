rule run_codeml:
    input:
        no_gap_file = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_nogap_msa.fasta",
        tree = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.bestTree_unrooted",
    params:
        wdir = "{workdir}/results/{query_fasta}/4_codeml",
    output:
        codeml_out="{workdir}/results/{query_fasta}/4_codeml/{query_fasta}_codeml",
    conda:
        "../envs/paml.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_run_codeml.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_run_codeml.out"
    cache: True
    resources:
        time=1200
    shell:
        "PATH={config[workdir]}/resources/paml4.9j/bin:$PATH && python scripts/codeml.py {input.no_gap_file} {input.tree} {output.codeml_out} {params.wdir} {config[seqtype]} {config[verbose]} {config[noisy]} {config[clock]} {config[aa_dist]} {config[aa_rate_file]} {config[model]} {config[icode]} {config[Mgene]} {config[fix_alpha]} {config[alpha]} {config[Malpha]} {config[ncatG]} {config[getSE]} {config[RateAncestor]} {config[Small_Diff]} {config[cleandata]} {config[method]} 2> {log}"
