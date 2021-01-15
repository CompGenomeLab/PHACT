rule compute_score:
    input:
        tree_file = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.bestTree_unrooted",
        codeml_out ="{output_folder}/results/{msa_name}/4_codeml/{msa_name}_codeml",
    output:
        "{output_folder}/results/{msa_name}/5_scores/{msa_name}_wol_param_{pattern}.csv",
        "{output_folder}/results/{msa_name}/5_scores/{msa_name}_wl_param_{pattern}.csv",
    params:
        out = "{output_folder}/results/{msa_name}/5_scores/{msa_name}",
        rst = "{output_folder}/results/{msa_name}/4_codeml/rst",
        fasta = "{output_folder}/resources/msa_files/{msa_name}.fasta",
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_{pattern}_compute_score.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_{pattern}_compute_score.out"
    conda:
        "../envs/r-base.yml"
    shell:
        "human_id=`python scripts/get_human_id.py {params.fasta}` && Rscript scripts/Blosum_compute_score.R {input.tree_file} {params.rst} {params.out} $human_id {config[weights]} 2>{log}"
