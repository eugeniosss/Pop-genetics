### PARSE TARGETS FILE BY CHR AND CREATE INTERMEDIARY TARGETS FILE PER CHR

rule parse:
    input:
        target=config["target_list"],
    output:
        temp("split_bed/{chr}.bed")
    log:
        "split_bed/{chr}.log"
    shell:
        """
        (awk '$1 == "{wildcards.chr}"' {input.target} > {output}) > {log} 2>&1
        """
