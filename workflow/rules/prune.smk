rule prune:
    input:
        blastp_out_file = "{workdir}/results/{query_id}/1_psiblast/{query_id}_blasthits.out",
        ml_tree = "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree",
        msa_file = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa.fasta",
        query_fasta = "{workdir}/results/{query_id}/1_psiblast/{query_id}.fasta",
    output:
        pruned_tree =  "{workdir}/results/{query_id}/6_pruned_msa_tree/{query_id}.pruned",
        pruned_msa = "{workdir}/results/{query_id}/6_pruned_msa_tree/{query_id}_pruned_nogap_msa.fasta",
    log:
        "{workdir}/workflow/logs/rules/{query_id}_pruned.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_pruned.out"
    cache: True
    conda:
        "../envs/prune.yml"
    shell:
        "python3 scripts/prune_msa_tree.py {input.blastp_out_file} {config[blast_hit_number]} {config[max_e_value]} {config[min_identity]} {config[max_identity]} {input.ml_tree} {output.pruned_tree} {input.msa_file} {output.pruned_msa} {input.query_fasta}"
