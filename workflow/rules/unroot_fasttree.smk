rule unroot_fasttree:
    input:
        input_tree = "{workdir}/results/{query_id}/6_fasttree/{query_id}.nwk",
        no_gap_msa_file = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa_fasttree.fasta",
    output:
      "{workdir}/results/{query_id}/6_fasttree/{query_id}.nwk_unrooted",
    log:
        "{workdir}/workflow/logs/rules/{query_id}_unroot_tree.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_unroot_tree.out"
    conda:
        "../envs/remove_outlier.yml"
    shell:
        "python  scripts/unroot_tree.py {input.input_tree}  {input.no_gap_msa_file} 2> {log}"
