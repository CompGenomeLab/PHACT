rule query_fasta:
    input:
        blastdb = expand("{workdir}/{blastdb_folder}/{blastdb}", workdir=config["workdir"], blastdb_folder=config["blastdb_folder"], blastdb=config["blastdb_file"])
    output:
        fasta_file= "{workdir}/results/{query_id}/1_psiblast/{query_id}.fasta"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_query_fasta.err" 
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_query_fasta.out"
    params:
        query_id = "{query_id}"
    cache: True
    conda:
        "../envs/python.yml"   
    shell:
        "python scripts/make_query_fasta.py {params.query_id}  {input.blastdb} {output.fasta_file} 2> {log}"
