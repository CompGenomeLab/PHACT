rule compute_score:
    input:
        tree_file = "{workdir}/results/{query_fasta}/3_mltree/{query_fasta}.raxml.bestTree_unrooted",
        codeml_out ="{workdir}/results/{query_fasta}/4_codeml/{query_fasta}_codeml",
    output:
        "{workdir}/results/{query_fasta}/5_scores/{query_fasta}_wol_param_{pattern}.csv",
        "{workdir}/results/{query_fasta}/5_scores/{query_fasta}_wl_param_{pattern}.csv",
    params:
        out = "{workdir}/results/{query_fasta}/5_scores/{query_fasta}",
        rst = "{workdir}/results/{query_fasta}/4_codeml/rst",
        fasta = "{workdir}/resources/query_fasta/{query_fasta}.fasta",
    log:
        "{workdir}/workflow/logs/rules/{query_fasta}_{pattern}_compute_score.err"
    benchmark:
        "{workdir}/workflow/logs/benchmarks/{query_fasta}_{pattern}_compute_score.out"
    cache: True
    conda:
        "../envs/r-base.yml"
    shell:
        "human_id=`python scripts/get_human_id.py {params.fasta}` && Rscript scripts/Blosum_compute_score.R {input.tree_file} {params.rst} {params.out} $human_id {config[weights]} 2>{log}"
