rule pruned_raxmlng_ancestral:
    input:
        pruned_nogap_msa = "{workdir}/results/{query_id}/6_pruned_msa_tree/{query_id}_pruned_nogap_msa.fasta",
        pruned_unrooted_tree = "{workdir}/results/{query_id}/6_pruned_msa_tree/{query_id}.pruned_unrooted",        
    output:
        pruned_ancestral_probabilities = "{workdir}/results/{query_id}/7_pruned_raxmlng_ancestral/pruned_{query_id}.raxml.ancestralProbs",
        pruned_ancestral_states = "{workdir}/results/{query_id}/7_pruned_raxmlng_ancestral/pruned_{query_id}.raxml.ancestralStates",
        pruned_ancestralTree = "{workdir}/results/{query_id}/7_pruned_raxmlng_ancestral/pruned_{query_id}.raxml.ancestralTree",

    params:
        pruned_raxml_ancestral_out_name = "{workdir}/results/{query_id}/7_pruned_raxmlng_ancestral/pruned_{query_id}",
    log:
        "{workdir}/workflow/logs/rules/{query_id}_pruned_raxmlng_ancestral.err"
    cache: True
    conda:
        "../envs/raxml-ng.yml"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_pruned_raxmlng_ancestral.out"
    shell:
        "raxml-ng --ancestral --msa {input.pruned_nogap_msa} --tree {input.pruned_unrooted_tree} --model {config[raxmlng_ancestral_model]} --prefix {params.pruned_raxml_ancestral_out_name}  --threads {config[raxml_ancestral_threads]}"


