rule unroot_pruned_tree:
    input:
        input_tree = "{workdir}/results/{query_id}/6_pruned_msa_tree/{query_id}.pruned",
    output:
        "{workdir}/results/{query_id}/6_pruned_msa_tree/{query_id}.pruned_unrooted",
    log:
        "{workdir}/workflow/logs/rules/{query_id}_pruned_unroot_tree.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_pruned_unroot_tree.out"
    conda:
        "../envs/r-base.yml"
    shell:
        "Rscript scripts/unroot_tree.R {input.input_tree} 2> {log}"
