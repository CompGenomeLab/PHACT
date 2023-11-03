rule iqtree_ancestral:
    input:
        no_gap_masked_msa = "{workdir}/results/{query_id}/2_msa/{query_id}_nogap_masked_msa.fasta",
        unrooted_tree = "{workdir}/results/{query_id}/6_fasttree/{query_id}.nwk_unrooted",        
    output:
        ancestral_probabilities = "{workdir}/results/{query_id}/7_iqtree_ancestral/{query_id}_nogap_msa.fasta.state",
        ancestralTree = "{workdir}/results/{query_id}/7_iqtree_ancestral/{query_id}_nogap_msa.fasta.treefile",
    params:
        iqtree_ancestral_out_name = "{workdir}/results/{query_id}/7_iqtree_ancestral/{query_id}_nogap_msa.fasta",
    conda:
        "../envs/iqtree.yml"
    log:
        "{workdir}/workflow/logs/rules/{query_id}_iqtree_ancestral.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_id}_iqtree_ancestral.out"
    resources:
        cpus=4
    shell:
        "iqtree2 -redo -s {input.no_gap_masked_msa} -te {input.unrooted_tree} -m {config[iqtree_ancestral_model]} -asr -safe -nt {resources.cpus} --prefix {params.iqtree_ancestral_out_name}"

