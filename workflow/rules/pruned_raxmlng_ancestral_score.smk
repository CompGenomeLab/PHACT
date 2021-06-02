rule pruned_raxml_anc_score:
    input:
        tree_file = "{workdir}/results/{query_id}/7_pruned_raxmlng_ancestral/pruned_{query_id}.raxml.ancestralTree",
        probabilities = "{workdir}/results/{query_id}/7_pruned_raxmlng_ancestral/pruned_{query_id}.raxml.ancestralProbs",
    output:
        "{workdir}/results/{query_id}/8_pruned_raxmlng_ancestral_scores/pruned_{query_id}_wol_param_{pattern}.csv",
        "{workdir}/results/{query_id}/8_pruned_raxmlng_ancestral_scores/pruned_{query_id}_wl_param_{pattern}.csv",
    params:
        out = "{workdir}/results/{query_id}/8_pruned_raxmlng_ancestral_scores/pruned_{query_id}",
        fasta = "{workdir}/results/{query_id}/6_pruned_msa_tree/{query_id}_pruned_nogap_msa.fasta",
        query_fasta = "{workdir}/results/{query_id}/1_psiblast/{query_id}.fasta",
    log:
        "{workdir}/workflow/logs/rules/{query_id}_pruned_raxmlanc_{pattern}_compute_score.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_pruned_raxmlanc_{pattern}_compute_score.out"
    cache: True
    conda:
        "../envs/r-base.yml"
    shell:
        "query=`python scripts/get_query.py {params.query_fasta}` && Rscript scripts/compute_score_RaxmlNg_Final.R {input.tree_file} {input.probabilities} {params.fasta} {params.out} $query {config[weights]} 2>{log}"


