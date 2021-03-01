rule get_blasthits:
    input:
        blastp_out = "{workdir}/results/{query_id}/1_psiblast/{query_id}_blasthits.out",
        query_fasta = "{workdir}/results/{query_id}/1_psiblast/{query_id}.fasta",
        blastdb = expand("{workdir}/{blastdb_folder}/{blastdb}", workdir=config["workdir"], blastdb_folder=config["blastdb_folder"], blastdb=config["blastdb_file"])
    output:
        blast_fasta = "{workdir}/results/{query_id}/1_psiblast/{query_id}_blasthits.fasta",
    resources:
        time_min=300
    conda:
        "../envs/python.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_get_blasthits.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_get_blasthits.out"
    cache: True
    shell:
        "python3 scripts/parse_blastp.py {input.blastp_out} {config[blast_hit_number]} {config[max_e_value]} {config[min_identity]} {input.blastdb} {input.query_fasta} 2> {log}"
