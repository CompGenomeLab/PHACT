rule get_blasthits:
    input:
        blastp_out = "{workdir}/results/{query_fasta}/1_blastp/{query_fasta}_blasthits.out",
	blastdb = expand("{workdir}/{blastdb_folder}/{blastdb}", workdir=config["workdir"], blastdb_folder=config["blastdb_folder"], blastdb=config["blastdb_file"]),
        query_fasta = "{workdir}/resources/msa_files/{query_fasta}.fasta",
    output:
        "{workdir}/results/{query_fasta}/1_blastp/{query_fasta}_blasthits.fasta",
    resources:
        time_min=300
    conda:
        "../envs/python.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_get_blasthits.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_get_blasthits.out"
    shell:
        "python3 scripts/parse_blastp.py {input.blastp_out} {config[blast_hit_number]} {config[max_e_value]} {config[min_identity]}  {config[max_identity]}  {input.blastdb} {input.query_fasta}  2> {log}"
