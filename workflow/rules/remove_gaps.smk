rule remove_gaps:
    input:
        query_file ="{workdir}/resources/msa_files/{query_fasta}.fasta",
        msa_file = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_blasthits_new_header_msa.fasta",
    output:
        no_gap_file = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_nogap_msa.fasta",
    conda:
        "../envs/python.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_remove_gaps.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_remove_gaps.out"
    shell:
        "python3 scripts/remove_gaps.py {input.query_file} {input.msa_file} 2> {log}"
