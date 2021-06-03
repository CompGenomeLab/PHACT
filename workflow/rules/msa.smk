rule msa:
    input:
        fasta = "{workdir}/results/{query_id}/1_psiblast/{query_id}_blasthits.fasta",
    output:
        msa_file = "{workdir}/results/{query_id}/2_msa/{query_id}_blasthits_msa.fasta",
    conda:
        "../envs/mafft.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_msa.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_msa.out"
    cache: True
    resources:
        cpus=4
    shell:
        "mafft{config[mafft_method]} --thread {resources.cpus} {input.fasta} > {output.msa_file} 2> {log}"
