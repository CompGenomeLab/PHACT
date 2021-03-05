rule remove_gaps:
    input:
        query_file = "{workdir}/results/{query_id}/1_psiblast/{query_id}.fasta",
        msa_file = "{workdir}/results/{query_id}/2_msa/{query_id}_blasthits_msa.fasta",
        tree_file = "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree",
    output:
        no_gap_file = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa.fasta",
    conda:
        "../envs/prune.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_remove_gaps.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_remove_gaps.out"
    cache: True
    shell:
        "python3 scripts/remove_gaps.py {input.query_file} {input.msa_file} {input.tree_file}  {output.no_gap_file} 2> {log}"
