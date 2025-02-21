def get_chromosomes_from_bed(bed_file):
    with open(bed_file) as f:
        return sorted(set(line.split()[0] for line in f if line.strip()))

# Use this function to extract chromosomes
chromosomes = get_chromosomes_from_bed(config["target_list"])

rule merge:
    input:
        # Dynamically create a list of VCF files for each chromosome (from the rule's outputs)
        expand("{{software}}/{chr}_variants.vcf", chr=chromosomes)  # Or use a list of chromosomes dynamically
    output:
        "{software}/variants.vcf"  # Output merged VCF file
    log:
        "{software}/final_merge.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    shell:
        """
        bcftools concat -o {output} {input} > {log} 2>&1
        """