rule psiblastp:
    input:
        blastdb = expand("{workdir}/{blastdb_folder}/{blastdb}", workdir=config["workdir"], blastdb_folder=config["blastdb_folder"], blastdb=config["blastdb_file"]),
        query_fasta = "{workdir}/results/{query_id}/1_psiblast/{query_id}.fasta"
    output:
        outfile = "{workdir}/results/{query_id}/1_psiblast/{query_id}_blasthits.out"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_psiblastp.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_psiblastp.out"
    cache: True
    conda:
        "../envs/blastp.yml"
    shell:
        "psiblast -query {input.query_fasta} -db {input.blastdb} -outfmt {config[outfmt]} -out {output.outfile} -max_target_seqs {config[max_target_seqs]} -num_iterations {config[num_iterations]}  2> {log}"
