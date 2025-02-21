### PARSE TARGETS FILE BY CHR AND CREATE INTERMEDIARY TARGETS FILE PER CHR

rule parse:
    input:
        target=config["target_list"],
    output:
        temp("split_bed/{chr}.bed")
    shell:
        """
        mkdir -p split_bed
        awk '$1 == "{wildcards.chr}"' {input.target} > {output}
        """
