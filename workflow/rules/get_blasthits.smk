rule get_blasthits:
    input:
        blastp_out = "{output_folder}/results/{msa_name}/1_blastp/{msa_name}_blasthits.out",
	blastdb = expand("{output_folder}/{blastdb_folder}/{blastdb}", output_folder=config["output_folder"], blastdb_folder=config["blastdb_folder"], blastdb=config["blastdb_file"]),
        query_fasta = "{output_folder}/resources/msa_files/{msa_name}.fasta",
    output:
        "{output_folder}/results/{msa_name}/1_blastp/{msa_name}_blasthits.fasta",
    resources:
        time_min=300
    conda:
        "../envs/python.yml"
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_get_blasthits.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_get_blasthits.out"
    shell:
        "python3 scripts/parse_blastp.py {input.blastp_out} {config[blast_hit_number]} {config[max_e_value]} {config[min_identity]}  {config[max_identity]}  {input.blastdb} {input.query_fasta}  2> {log}"
