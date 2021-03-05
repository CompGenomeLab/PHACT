rule raxmlng_ancestral:
    input:
        nogap_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_msa.fasta",
        unrooted_tree = "{workdir}/results/{query_id}/3_mltree/{query_id}.raxml.bestTree_unrooted",        
    output:
        ancestral_probabilities = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralProbs",
        ancestral_states = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralStates",
        ancestralTree = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}.raxml.ancestralTree",
    params:
        raxml_ancestral_out_name = "{workdir}/results/{query_id}/4_raxmlng_ancestral/{query_id}",
    conda:
        "../envs/raxml-ng.yml"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_raxmlng_ancestral.out"
    cache: True
    resources:
        time_min=7200,cpus=2
    shell:
        "raxml-ng --ancestral --msa {input.nogap_msa} --tree {input.unrooted_tree} --model {config[raxmlng_ancestral_model]} --prefix {params.raxml_ancestral_out_name} --threads {resources.cpus}"
