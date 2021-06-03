rule raxmlng_ancestral:
    input:
        no_gap_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa.fasta",
        unrooted_tree = "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree_unrooted",        
    output:
        ancestral_probabilities = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralProbs",
        ancestral_states = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralStates",
        ancestralTree = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralTree",
    params:
        raxml_ancestral_out_name = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}",
    conda:
        "../envs/raxml-ng.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_raxmlng_ancestral.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_raxmlng_ancestral.out"
    cache: True
    resources:
        cpus=4
    shell:
        "raxml-ng --ancestral --msa {input.no_gap_msa} --tree {input.unrooted_tree} --model {config[raxmlng_ancestral_model]} --prefix {params.raxml_ancestral_out_name} --threads {resources.cpus}"
