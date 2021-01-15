rule msa:
    input:
        fasta = "{output_folder}/results/{msa_name}/1_blastp/{msa_name}_blasthits_new_header.fasta",
    output:
        msa_file = "{output_folder}/results/{msa_name}/2_msa/{msa_name}_blasthits_new_header_msa.fasta",
    conda:
        "../envs/mafft.yml"
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_msa.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_msa.out"
    shell:
        "mafft{config[mafft_method]} {input.fasta} > {output.msa_file} 2> {log}"
