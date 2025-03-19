
rule create_keep_file:
    input:
         "{subset}.txt"
    output:
         temp("{subset}_mod.txt")
    shell:
         "awk '{{print $1,$1,0,0,0,-9}}' {input} > {output}"

rule vcf_to_bed_with_plink:
    input:
        VCF=config["initial_vcf"],
        KEEP="{subset}_mod.txt"
    output:
        BED=temp("plink_{subset}/{subset}.bed"),
        BIM=temp("plink_{subset}/{subset}.bim"),
        FAM=temp("plink_{subset}/{subset}.fam"),
        FRQ=temp("plink_{subset}/{subset}.frq"),
        MAP=temp("plink_{subset}/{subset}.map"),
        NOSEX=temp("plink_{subset}/{subset}.nosex"),
        PED=temp("plink_{subset}/{subset}.ped"),
    log:
        "plink_{subset}/{subset}_snake.log"
    threads: 1
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        prefix="plink_{subset}/{subset}",
    shell:
        "(plink \
            --vcf {input.VCF} \
            --biallelic-only strict \
	        --keep {input.KEEP}\
            --snps-only \
            --set-missing-var-ids @_# \
            --make-bed \
            --allow-extra-chr \
            --double-id \
            --freq \
            --recode \
            --out {params.prefix}) > {log} 2>&1"