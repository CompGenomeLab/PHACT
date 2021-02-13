rule header_update:
    input:
        blasthit_file = "{workdir}/results/{query_fasta}/1_blastp/{query_fasta}_blasthits.fasta",
    output:
        new_header_fasta = "{workdir}/results/{query_fasta}/1_blastp/{query_fasta}_blasthits_new_header.fasta",
    conda:
        "../envs/python.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_header_update.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_header_update.out"
    cache: True
    shell:
        "python3 scripts/change_fasta_header.py {input.blasthit_file} {output.new_header_fasta} 2> {log}"
