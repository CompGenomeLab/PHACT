rule msa_masking_fasttree:
    input:
        no_gap_file = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa_fasttree.fasta",
    output:
        no_gap_masked_file = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_masked_msa.fasta",
    params:
        query_id = "{query_id}"
    conda:
        "../envs/r-base.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_masked_msa.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_masked_msa.out"
    cache: True
    shell:
        "Rscript scripts/MSA_Masking.R {input.no_gap_file} {params.query_id} {output.no_gap_masked_file} 2> {log}"
