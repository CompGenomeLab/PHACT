rule fast_tree:
    input:
        trimmed_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_trimmed_msa.fasta",        
    output:
        bestTree = "{workdir}/results/{query_id}/6_fasttree/{query_id}.nwk"
    params:
        raxml_out_name = "{query_id}",
    conda:
        "../envs/fasttree.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_fastreetree.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_fasttree_tree.out"
    cache: True
    shell:
        "FastTreeMP {config[fasttree_model]} {input.trimmed_msa}  >  {output.bestTree} 2>{log}"

