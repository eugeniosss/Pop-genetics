### PARSE TARGETS FILE BY CHR AND CREATE INTERMEDIARY TARGETS FILE PER CHR

rule list_chromosomes:
    input:
        config["target_list"],
    output:
        temp("chromosomes.txt")
    shell:
        """
        cut -f1 {input} | uniq > {output}
        """

rule parse:
    input:
        target=config["target_list"],
        list="chromosomes.txt"
    output:
        temp("split_bed/{chr}.bed")
    shell:
        """
        mkdir -p split_bed
        awk '$1 == "{wildcards.chr}"' {input.target} > {output}
        """
