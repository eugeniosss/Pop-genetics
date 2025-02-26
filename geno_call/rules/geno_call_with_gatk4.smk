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

def get_bams_dict(wildcards):
    return create_bams_dict(f"{wildcards.dataset}.txt")

def get_chromosomes_from_bed(bed_file):
    with open(bed_file) as f:
        return sorted(set(line.split()[0] for line in f if line.strip()))

#bams_dict=create_bams_dict(config["bams"])
#bams_dict = lambda wildcards: create_bams_dict(f"{wildcards.dataset}.txt")
bams_dict = lambda wildcards: create_bams_dict(f"{wildcards.dataset}.txt")

chromosomes = get_chromosomes_from_bed(config["target_list"])


rule haplo_caller_with_gatk:
    input:
        reference=config["ref_genome"],
        bam=lambda wildcards: bams_dict(wildcards)[wildcards.SAMPLE],
        bed_file="split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        gvcf=temp("gatk_{dataset}/{chr}_{SAMPLE}_variants.gvcf"),
        idx=temp("gatk_{dataset}/{chr}_{SAMPLE}_variants.gvcf.idx")
    log:
        "gatk_{dataset}/{chr}_{SAMPLE}.log"
    threads: 1
    conda:
        config["dir"] + "envs/gatk.yml"
    params:
        options=config["options"]["gatk_HaplotypeCaller"],
        mem_mb=int((config["options"]["gatk_HaplotypeCaller"]).split("--java-options -Xmx")[1].split("m")[0])*1.2
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
            {params.options} \
            -O {output.gvcf}) > {log} 2>&1
        """

rule CombineGVCFs_with_gatk:
    input:
        reference=config["ref_genome"],
        gvcfs=lambda wildcards: expand("gatk_{dataset}/{chr}_{SAMPLE}_variants.gvcf", 
            dataset=wildcards.dataset, 
            chr=wildcards.chr,
            SAMPLE=list(get_bams_dict(wildcards).keys()))
    output:
        gvcf=temp("gatk_{dataset}/{chr}_variants.gvcf"),
        idx=temp("gatk_{dataset}/{chr}_variants.gvcf.idx")
    log:
        "gatk_{dataset}/{chr}_merging_GVCFs.log"
    params:
        variants=lambda wildcards: create_variant_string(expand("gatk_{dataset}/{chr}_{SAMPLE}_variants.gvcf", 
           dataset=wildcards.dataset, 
           chr=wildcards.chr, 
           SAMPLE=list(get_bams_dict(wildcards).keys()))),
        options=config["options"]["gatk_CombineGVCFs"],
        mem_mb=int((config["options"]["gatk_CombineGVCFs"]).split("--java-options -Xmx")[1].split("m")[0])*1.2
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
            {params.options} \
            -O {output.gvcf}) > {log} 2>&1
        """

rule GenotypeGVCFs_with_gatk:
    input:
        reference=config["ref_genome"],
        gvcf="gatk_{dataset}/{chr}_variants.gvcf"
    output:
        vcf=temp("gatk_{dataset}/{chr}_variants.vcf"),
        idx=temp("gatk_{dataset}/{chr}_variants.vcf.idx")
    log:
        "gatk_{dataset}/{chr}_genotyping.log"
    threads: 1
    conda:
        config["dir"] + "envs/gatk.yml"
    params:
        options=config["options"]["gatk_GenotypeGVCFs"],
        mem_mb=int((config["options"]["gatk_GenotypeGVCFs"]).split("--java-options -Xmx")[1].split("m")[0])*1.2
    wildcard_constraints:
        chr="|".join(chromosomes)  # Restrains 'chr' to values from the file
    shell:
        """
        (gatk \
            GenotypeGVCFs \
            -R {input.reference} \
            -V {input.gvcf}\
            {params.options} \
            -O {output.vcf}) > {log} 2>&1
        """