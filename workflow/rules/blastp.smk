rule blastp:
    input:
        fasta = "{output_folder}/resources/msa_files/{msa_name}.fasta",
        blastdb = expand("{output_folder}/{blastdb_folder}/{blastdb}", output_folder=config["output_folder"], blastdb_folder=config["blastdb_folder"], blastdb=config["blastdb_file"])
    output:
        outfile = "{output_folder}/results/{msa_name}/1_blastp/{msa_name}_blasthits.out"
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_blastp.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_blastp.out"
    conda:
        "../envs/blastp.yml"
    shell:
        "blastp -query {input.fasta} -db {input.blastdb} -outfmt {config[outfmt]} -out {output.outfile} -max_target_seqs {config[max_target_seqs]} 2> {log}"
