
rule filter_with_plink:
    input:
        BED="plink_{subset}/{subset}_no_trans_pseudo_LD.bed",
        BIM="plink_{subset}/{subset}_no_trans_pseudo_LD.bim",
        FAM="plink_{subset}/{subset}_no_trans_pseudo_LD.fam",
        NOSEX="plink_{subset}/{subset}_no_trans_pseudo_LD.nosex",
    output:
        BED="admixture_{subset}/{subset}.bed",
        BIM="admixture_{subset}/{subset}.bim",
        FAM="admixture_{subset}/{subset}.fam",
    log:
        "admixture_{subset}/{subset}_snake.log"
    threads: 1
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        input_prefix="plink_{subset}/{subset}_no_trans_pseudo_LD",
        output_prefix="admixture_{subset}/{subset}"
    shell:
        """
        (plink \
        --bfile {params.input_prefix} \
        --allow-extra-chr \
        --geno 0.3 \
        --make-bed \
        --out {params.output_prefix}

        awk '{{$1=1 ; print ;}}' {params.output_prefix}.bim > {params.output_prefix}_tmp.bim

        mv {params.output_prefix}_tmp.bim {params.output_prefix}.bim

        ) > {log} 2>&1
        """

rule run_admixture:
    input:
        "admixture_{subset}/{subset}.bed",
    output:
        P="admixture_{subset}/{subset}.{K}.P",
        Q="admixture_{subset}/{subset}.{K}.Q"
    log:
        "admixture_{subset}/{subset}.{K}.log"
    conda:
        config["dir"] + "envs/admix.yml"
    params:
        K="{K}",
        subset="{subset}",
    threads: 1
    shell:
        """
        (cd admixture_{params.subset}
        admixture --cv -j1 --haploid="*" {params.subset}.bed {params.K}) > {log} 2>&1
        """