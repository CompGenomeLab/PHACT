rule header_update:
    input:
        blasthit_file = "{output_folder}/results/{msa_name}/1_blastp/{msa_name}_blasthits.fasta",
    output:
        new_header_fasta = "{output_folder}/results/{msa_name}/1_blastp/{msa_name}_blasthits_new_header.fasta",
    conda:
        "../envs/python.yml"
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_header_update.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_header_update.out"
    shell:
        "python3 scripts/change_fasta_header.py {input.blasthit_file} {output.new_header_fasta} 2> {log}"
