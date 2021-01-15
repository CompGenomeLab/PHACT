rule ml_tree:
    input:
        trimmed_msa = "{output_folder}/results/{msa_name}/2_msa/{msa_name}_trimmed_msa.fasta",        
    output:
        bestTree = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.bestTree",
        log_file = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.log",
        bestModel = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.bestModel",
        mlTrees = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.mlTrees",
        rba = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.rba",
        startTree = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}.raxml.startTree",
    params:
        raxml_out_name = "{output_folder}/results/{msa_name}/3_mltree/{msa_name}",
    conda:
        "../envs/raxml-ng.yml"
    benchmark:
        "{output_folder}/workflow/logs/benchmarks/{msa_name}_ml_tree.out"
    resources:
        time_min=7200,cpus=12
    shell:
        "raxml-ng  --msa {input.trimmed_msa} --data-type AA --model {config[raxml_model]} --prefix {params.raxml_out_name} --seed {config[raxml_seed]} --threads {config[raxml_threads]}"
