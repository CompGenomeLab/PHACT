rule remove_outliers:
    input:
        nogap_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa.fasta",
        input_tree = "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree"        
    output:
        no_outlier_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_no_outlier_no_gap.fasta",
        no_outlier_tree = "{workdir}/results/{query_id}/3_mltree/{query_id}.no_outlier.nwk",
  
    conda:
        "../envs/remove_outlier.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_remove_outliers.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_no_outlier.out"
    resources:
        time_min=7200,cpus=1
    shell:
        "python3 scripts/remove_outliers.py {input.input_tree} {input.nogap_msa} {output.no_outlier_tree} {output.no_outlier_msa} {config[max_deviations]} 2> {log}"

