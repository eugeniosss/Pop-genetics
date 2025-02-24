
rule create_keep_file:
    input:
         config["subsets"]+".txt"
    output:
         config["subsets"]+"_mod.txt"
    shell:
         "awk '{{print $1,$1,0,0,0,-9}}' {input} > {output}"

rule vcf_to_bed_with_plink:
    input:
        VCF=config["initial_vcf"],
        KEEP=config["subsets"]+"_mod.txt"
    output:
        BED=config["subsets"]+".bed",
        BIM=config["subsets"]+".bim",
        FAM=config["subsets"]+".fam",
        FRQ=config["subsets"]+".frq",
        MAP=config["subsets"]+".map",
        NOSEX=config["subsets"]+".nosex",
        PED=config["subsets"]+".ped",
    log:
        config["subsets"]+"snake.log"
    threads: 1
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        prefix=config["subsets"],
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