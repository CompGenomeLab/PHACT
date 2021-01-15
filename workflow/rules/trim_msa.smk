rule trim_msa:
    input:
        msa_file = "{output_folder}/results/{msa_name}/2_msa/{msa_name}_blasthits_new_header_msa.fasta",
    output:
        trimmed_msa = "{output_folder}/results/{msa_name}/2_msa/{msa_name}_trimmed_msa.fasta",
    conda:
        "../envs/trimal.yml"
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_trim_msa.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_trim_msa.out"
    shell:
        "trimal -in {input.msa_file} -out {output.trimmed_msa} {config[trimal_method]} -keepseqs 2> {log}"
