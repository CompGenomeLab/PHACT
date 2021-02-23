rule ml_tree:
    input:
        trimmed_msa = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_trimmed_msa.fasta",        
    output:
        bestTree = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.bestTree",
    params:
        raxml_out_name = "{query_fasta}",
    conda:
        "../envs/raxml-ng.yml"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_ml_tree.out"
    cache: True
    resources:
        time_min=7200,cpus=16
    shell:
        "scripts/raxml-ng.sh {input.trimmed_msa}  {config[raxml_model]} {params.raxml_out_name} {config[raxml_seed]} {output.bestTree} {config[raxml_tree_number]}"
