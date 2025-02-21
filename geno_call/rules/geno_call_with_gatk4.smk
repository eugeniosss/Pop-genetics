def create_variant_string(input_list):
    # Prefix each entry in the list with '--variant' and join them into a single string
    return " ".join(f"--variant {item}" for item in input_list)

def create_bams_dict(bams_file):
    with open(bams_file) as f:
        bams_dict = {
            os.path.splitext(os.path.basename(line.strip()))[0]: line.strip()
            for line in f
        }
    return bams_dict

def get_chromosomes_from_bed(bed_file):
    with open(bed_file) as f:
        return sorted(set(line.split()[0] for line in f if line.strip()))

bams_dict=create_bams_dict(config["bams"])

chromosomes = get_chromosomes_from_bed(config["target_list"])


rule haplo_caller_with_gatk:
    input:
        reference=config["ref_genome"],
        bam=lambda wildcards: bams_dict[wildcards.SAMPLE],
        bed_file="split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        "gatk/{chr}_{SAMPLE}_variants.gvcf"   # Output VCF per chromosome
    log:
        "gatk/{chr}_{SAMPLE}.log"
    threads: 1
    conda:
        config["dir"] + "envs/gatk.yml"
    wildcard_constraints:
        chr="|".join(chromosomes)  # Restrains 'chr' to values from the file
    shell:
        """
        (gatk \
            HaplotypeCaller \
            -R {input.reference} \
            --native-pair-hmm-threads {threads} \
            -I {input.bam} \
            -L {input.bed_file} \
            -ERC GVCF \
            -O {output}) > {log} 2>&1
        """

rule CombineGVCFs_with_gatk:
    input:
        reference=config["ref_genome"],
#        gvcfs = ["gatk/{{chr}}_{SAMPLE}_variants.gvcf".format(SAMPLE=SAMPLE) for SAMPLE in list(bams_dict.keys())]
        gvcfs = expand("gatk/{{chr}}_{SAMPLE}_variants.gvcf", SAMPLE=list(bams_dict.keys()))  # Or use a list of chromosomes dynamically
    output:
        "gatk/{chr}_variants.gvcf"
    log:
        "gatk/{chr}_merging_GVCFs.log"
    params:
        variants=create_variant_string(expand("gatk/{{chr}}_{SAMPLE}_variants.gvcf", SAMPLE=list(bams_dict.keys()))),
    threads: 1
    conda:
        config["dir"] + "envs/gatk.yml"
    wildcard_constraints:
        chr="|".join(chromosomes)  # Restrains 'chr' to values from the file
    shell:
        """
        (gatk \
            CombineGVCFs \
            -R {input.reference} \
            {params.variants} \
            -O {output}) > {log} 2>&1
        """

rule GenotypeGVCFs_with_gatk:
    input:
        reference=config["ref_genome"],
        gvcf="gatk/{chr}_variants.gvcf"
    output:
        "gatk/{chr}_variants.vcf"
    log:
        "gatk/{chr}_genotyping.log"
    threads: 1
    conda:
        config["dir"] + "envs/gatk.yml"
    shell:
        """
        (gatk \
            GenotypeGVCFs \
            -R {input.reference} \
            -V {input.gvcf}\
            -O {output}) > {log} 2>&1
        """
