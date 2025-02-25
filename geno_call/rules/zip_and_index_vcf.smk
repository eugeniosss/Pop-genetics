rule zip_vcf:
    input:
        "{software}_{dataset}/variants.vcf"  # Input from final_merge rule
    output:
        "{software}_{dataset}/variants.vcf.gz"  # Output gzipped VCF file
    log:
        "{software}_{dataset}/zip_vcf.log"
    threads: 1
    shell:
        """
        bgzip -c {input} > {output} 2> {log}
        """

rule index_vcf_gz:
    input:
        "{software}_{dataset}/variants.vcf.gz"  # Input gzipped VCF file
    output:
        "{software}_{dataset}/variants.vcf.gz.tbi"  # Output index file
    threads: 1
    log:
        "{software}_{dataset}/index_vcf.log"
    shell:
        """
        tabix -p vcf {input} > {log} 2>&1
        """
