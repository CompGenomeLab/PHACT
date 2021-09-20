rule unroot_fasttree:
    input:
        input_tree = "{workdir}/results/{query_id}/6_fasttree/{query_id}.nwk"
    output:
      "{workdir}/results/{query_id}/6_fasttree/{query_id}.nwk_unrooted",
    log:
        "{workdir}/workflow/logs/rules/{query_id}_unroot_tree.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_unroot_tree.out"
    conda:
        "../envs/r-base.yml"
    shell:
        "Rscript scripts/unroot_tree.R {input.input_tree} 2> {log}"

