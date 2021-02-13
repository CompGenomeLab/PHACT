rule msa:
    input:
        fasta = "{workdir}/results/{query_fasta}/1_blastp/{query_fasta}_blasthits_new_header.fasta",
    output:
        msa_file = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_blasthits_new_header_msa.fasta",
    conda:
        "../envs/mafft.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_msa.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_msa.out"
    cache: True
    resources:
        cpus=8
    shell:
        "mafft{config[mafft_method]} --thread {resources.cpus} {input.fasta} > {output.msa_file} 2> {log}"
