
rule create_keep_file:
    input:
         config["subsets"]+".txt"
    output:
         temp(config["subsets"]+"_mod.txt")
    shell:
         "awk '{{print $1,$1,0,0,0,-9}}' {input} > {output}"

rule vcf_to_bed_with_plink:
    input:
        VCF=config["initial_vcf"],
        KEEP=config["subsets"]+"_mod.txt"
    output:
        BED=temp("plink/"+config["subsets"]+".bed"),
        BIM=temp("plink/"+config["subsets"]+".bim"),
        FAM=temp("plink/"+config["subsets"]+".fam"),
        FRQ=temp("plink/"+config["subsets"]+".frq"),
        MAP=temp("plink/"+config["subsets"]+".map"),
        NOSEX=temp("plink/"+config["subsets"]+".nosex"),
        PED=temp("plink/"+config["subsets"]+".ped"),
    log:
        "plink/"+config["subsets"]+"snake.log"
    threads: 1
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        prefix="plink/"+config["subsets"],
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