rule trim_msa:
    input:
        msa_file = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_blasthits_new_header_msa.fasta",
    output:
        trimmed_msa = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_trimmed_msa.fasta",
    conda:
        "../envs/trimal.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_trim_msa.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_trim_msa.out"
    cache: True
    shell:
        "trimal -in {input.msa_file} -out {output.trimmed_msa} {config[trimal_method]} -keepseqs 2> {log}"
