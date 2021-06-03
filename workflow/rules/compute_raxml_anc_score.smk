rule raxml_anc_score:
    input:
        tree_file = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralTree",
        probabilities = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralProbs",
    output:
        "{workdir}/results/{query_id}/5_raxmlng_ancestral_scores/{query_id}_wol_param_{pattern}.csv",
        "{workdir}/results/{query_id}/5_raxmlng_ancestral_scores/{query_id}_wl_param_{pattern}.csv",
    params:
        out = "{workdir}/results/{query_id}/5_raxmlng_ancestral_scores/{query_id}",
        fasta = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa.fasta",
        query_fasta = "{workdir}/results/{query_id}/1_psiblast/{query_id}.fasta" 
    log:
        "{workdir}/workflow/logs/rules/{query_id}_raxmlanc_{pattern}_compute_score.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_raxmlanc_{pattern}_compute_score.out"
    cache: True
    conda:
        "../envs/r-base.yml"
    shell:
        "query=`python scripts/get_query.py {params.query_fasta}` && Rscript scripts/compute_score_RaxmlNg_Final.R {input.tree_file} {input.probabilities} {params.fasta} {params.out} $query {config[weights]} 2>{log}"
