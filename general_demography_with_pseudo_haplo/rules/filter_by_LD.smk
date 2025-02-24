rule LD_filter_list:
    input:
        BED=config["subsets"]+"_no_trans_pseudo.bed",
        BIM=config["subsets"]+"_no_trans_pseudo.bim",
        FAM=config["subsets"]+"_no_trans_pseudo.fam",
    output:
        IN=config["subsets"]+"_no_trans_pseudo.prune.in",
        OUT=config["subsets"]+"_no_trans_pseudo.prune.out",
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        input_prefix=config["subsets"]+"_no_trans_pseudo",
        output_prefix=config["subsets"]+"_no_trans_pseudo",
    log:
        config["subsets"]+"_no_trans_pseudo_LD_snake.log"
    shell:
        "(plink \
        --bfile {params.input_prefix} \
        --indep-pairwise 200 25 0.6 \
        --allow-extra-chr \
        --out {params.output_prefix})  > {log} 2>&1"

rule LD_filter_with_plink:
    input:
        BED=config["subsets"]+"_no_trans_pseudo.bed",
        BIM=config["subsets"]+"_no_trans_pseudo.bim",
        FAM=config["subsets"]+"_no_trans_pseudo.fam",
        IN=config["subsets"]+"_no_trans_pseudo.prune.in",
    output:
        BED=config["subsets"]+"_no_trans_pseudo_LD.bed",
        BIM=config["subsets"]+"_no_trans_pseudo_LD.bim",
        FAM=config["subsets"]+"_no_trans_pseudo_LD.fam",
        NOSEX=config["subsets"]+"_no_trans_pseudo_LD.nosex",
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        input_prefix=config["subsets"]+"_no_trans_pseudo",
        output_prefix=config["subsets"]+"_no_trans_pseudo_LD",
    log:
        config["subsets"]+"_no_trans_pseudo_LD_snake.log"
    shell:
        "(plink \
        --bfile {params.input_prefix} \
        --extract {input.IN} \
        --make-bed \
        --allow-extra-chr \
        --out {params.output_prefix}) > {log} 2>&1"