
rule geno_call_with_freebayes:
    input:
        reference=config["ref_genome"],  # Reference genome
        bed_file="split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        "freebayes/{chr}_variants.vcf"   # Output VCF per chromosome
    log:
        "freebayes/{chr}.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    params:
        bams=config["bams"],  # BAM file list (define in config or Snakefile)
    shell:
        """
        mkdir -p freebayes
        freebayes \
        --fasta {input.reference} \
        -L {params.bams} \
        -t {input.bed_file} > {output} 2> {log}
        """
