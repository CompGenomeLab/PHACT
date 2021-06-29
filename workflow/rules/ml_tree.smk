rule ml_tree:
    input:
        trimmed_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_trimmed_msa.fasta",        
    output:
        bestTree = "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree"
    params:
        raxml_out_name = "{query_id}",
    conda:
        "../envs/raxml-ng.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_mltree.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_ml_tree.out"
    cache: True
    resources:
        cpus=8
    shell:
        "scripts/raxml-ng.sh {input.trimmed_msa}  {config[raxml_model]} {params.raxml_out_name} {config[raxml_seed]} {output.bestTree} {config[raxml_tree_number]} 2>{log}"
