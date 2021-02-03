rule unroot_tree:
    input:
        input_tree = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.bestTree",
    output:
      "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.bestTree_unrooted",
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_unroot_tree.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_unroot_tree.out"
    conda:
        "../envs/r-base.yml"
    shell:
        "Rscript scripts/unroot_tree.R {input.input_tree} 2> {log}"
