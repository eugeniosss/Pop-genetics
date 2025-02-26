def get_chromosomes_from_bed(bed_file):
    seen = set()
    chromosomes = []
    with open(bed_file) as f:
        for line in f:
            chrom = line.split()[0]
            if chrom not in seen and chrom.strip():
                chromosomes.append(chrom)
                seen.add(chrom)
    return chromosomes

# Use this function to extract chromosomes
chromosomes = get_chromosomes_from_bed(config["target_list"])

rule merge:
    input:
        # Dynamically create a list of VCF files for each chromosome (from the rule's outputs)
        expand("{{software}}_{{dataset}}/{chr}_variants.vcf", chr=chromosomes)  # Or use a list of chromosomes dynamically
        #lambda wildcards: [f"{{software}}/{chrom}_variants.vcf" for chrom in chromosomes]
    output:
        temp("{software}_{dataset}/variants_{software}_{dataset}.vcf")  # Output merged VCF file
    log:
        "{software}_{dataset}/final_merge.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    shell:
        """
        bcftools concat -o {output} {input} > {log} 2>&1
        echo {input} >> {log}
        """