rule trim_msa:
    input:
        msa_file = "{workdir}/results/{query_id}/2_msa/{query_id}_blasthits_msa.fasta",
    output:
        trimmed_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_trimmed_msa.fasta",
    conda:
        "../envs/trimal.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_trim_msa.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_trim_msa.out"
    cache: True
    shell:
        "trimal -in {input.msa_file} -out {output.trimmed_msa} {config[trimal_method]} 2> {log}"
