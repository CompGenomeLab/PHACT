rule unroot_tree:
    input:
        input_tree = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.bestTree",
    output:
      "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.bestTree_unrooted",
    log:
        "{output_folder}/workflow/logs/rules/{msa_name}_unroot_tree.err"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_unroot_tree.out"
    conda:
        "../envs/r-base.yml"
    shell:
        "Rscript scripts/unroot_tree.R {input.input_tree} 2> {log}"
