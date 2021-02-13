rule blastp:
    input:
        fasta = "{workdir}/resources/query_fasta/{query_fasta}.fasta",
        blastdb = expand("{workdir}/{blastdb_folder}/{blastdb}", workdir=config["workdir"], blastdb_folder=config["blastdb_folder"], blastdb=config["blastdb_file"])
    output:
        outfile = "{workdir}/results/{query_fasta}/1_blastp/{query_fasta}_blasthits.out"
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_blastp.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_blastp.out"
    cache: True
    conda:
        "../envs/blastp.yml"
    shell:
        "blastp -query {input.fasta} -db {input.blastdb} -outfmt {config[outfmt]} -out {output.outfile} -max_target_seqs {config[max_target_seqs]} 2> {log}"
