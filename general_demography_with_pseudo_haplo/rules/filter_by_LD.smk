rule LD_filter_list:
    input:
        BED="plink/"+config["subsets"]+"_no_trans_pseudo.bed",
        BIM="plink/"+config["subsets"]+"_no_trans_pseudo.bim",
        FAM="plink/"+config["subsets"]+"_no_trans_pseudo.fam",
    output:
        IN=temp("plink/"+config["subsets"]+"_no_trans_pseudo.prune.in"),
        OUT=temp("plink/"+config["subsets"]+"_no_trans_pseudo.prune.out"),
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        input_prefix="plink/"+config["subsets"]+"_no_trans_pseudo",
        output_prefix="plink/"+config["subsets"]+"_no_trans_pseudo",
    log:
        "plink/"+config["subsets"]+"_no_trans_pseudo_LD_snake.log"
    shell:
        "(plink \
        --bfile {params.input_prefix} \
        --indep-pairwise 200 25 0.6 \
        --allow-extra-chr \
        --out {params.output_prefix})  > {log} 2>&1"

rule LD_filter_with_plink:
    input:
        BED="plink/"+config["subsets"]+"_no_trans_pseudo.bed",
        BIM="plink/"+config["subsets"]+"_no_trans_pseudo.bim",
        FAM="plink/"+config["subsets"]+"_no_trans_pseudo.fam",
        IN="plink/"+config["subsets"]+"_no_trans_pseudo.prune.in",
    output:
        BED="plink/"+config["subsets"]+"_no_trans_pseudo_LD.bed",
        BIM="plink/"+config["subsets"]+"_no_trans_pseudo_LD.bim",
        FAM="plink/"+config["subsets"]+"_no_trans_pseudo_LD.fam",
        NOSEX="plink/"+config["subsets"]+"_no_trans_pseudo_LD.nosex",
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        input_prefix="plink/"+config["subsets"]+"_no_trans_pseudo",
        output_prefix="plink/"+config["subsets"]+"_no_trans_pseudo_LD",
    log:
        "plink/"+config["subsets"]+"_no_trans_pseudo_LD_snake.log"
    shell:
        "(plink \
        --bfile {params.input_prefix} \
        --extract {input.IN} \
        --make-bed \
        --allow-extra-chr \
        --out {params.output_prefix}) > {log} 2>&1"