rule unroot_tree:
    input:
        input_tree = "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree"
    output:
      "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree_unrooted",
    log:
        "{workdir}/workflow/logs/rules/{query_id}_unroot_tree.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_unroot_tree.out"
    conda:
        "../envs/r-base.yml"
    cache: True
    shell:
        "Rscript scripts/unroot_tree.R {input.input_tree} 2> {log}"
