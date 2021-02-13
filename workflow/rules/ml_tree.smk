rule ml_tree:
    input:
        trimmed_msa = "{workdir}/results/{query_fasta}/2_msa/{query_fasta}_trimmed_msa.fasta",        
    output:
        bestTree = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.bestTree",
        log_file = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.log",
        bestModel = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.bestModel",
        mlTrees = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.mlTrees",
        rba = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.rba",
        startTree = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.startTree",
    params:
        raxml_out_name = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}",
    conda:
        "../envs/raxml-ng.yml"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_ml_tree.out"
    cache: True
    resources:
        time_min=7200,cpus=12
    shell:
        "raxml-ng  --msa {input.trimmed_msa} --data-type AA --model {config[raxml_model]} --prefix {params.raxml_out_name} --seed {config[raxml_seed]} --threads {config[raxml_threads]}"
