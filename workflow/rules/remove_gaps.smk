rule remove_gaps:
    input:
        query_file ="{output_folder}/resources/msa_files/{msa_name}.fasta",
        msa_file = "{output_folder}/results/{msa_name}/2_msa/{msa_name}_blasthits_new_header_msa.fasta",
    output:
        no_gap_file = "{output_folder}/results/{msa_name}/2_msa/{msa_name}_nogap_msa.fasta",
    conda:
        "../envs/python.yml"
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_remove_gaps.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_remove_gaps.out"
    shell:
        "python3 scripts/remove_gaps.py {input.query_file} {input.msa_file} 2> {log}"
