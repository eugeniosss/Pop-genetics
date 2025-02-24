rule remove_trans:
    input:
        FRQ=config["subsets"]+".frq",
    output:
        FRQ=config["subsets"]+".snps",
    run:
        import csv

        IN_FREQ_FILE = input[0]
        OUT_FILE = output[0]
        transversion = [{'A', 'C'}, {'A', 'T'}, {'G', 'C'}, {'G', 'T'}]

        with open(OUT_FILE, 'w') as out_file:
            out_no_trans = csv.writer(out_file, delimiter='\t')

            with open(IN_FREQ_FILE) as tsv:
                for line in csv.reader(tsv, dialect="excel", delimiter=' ', skipinitialspace=True):
                    alleles = set(line[2] + line[3])
                    if alleles in transversion:
                        out_no_trans.writerow(line)

rule remove_trans_with_plink:
    input:
        BED=config["subsets"]+".bed",
        BIM=config["subsets"]+".bim",
        FAM=config["subsets"]+".fam",
        SNPS=config["subsets"]+".snps",
    output:
        MAP=config["subsets"]+"_no_trans.map",
        NOSEX=config["subsets"]+"_no_trans.nosex",
        PED=config["subsets"]+"_no_trans.ped",
        LOG=config["subsets"]+"_no_trans.log"
    log:
        config["subsets"]+"_no_trans_snake.log"
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        input_prefix=config["subsets"],
        output_prefix=config["subsets"]+"_no_trans",
    shell:
        "(plink \
            --bfile {params.input_prefix} \
            --extract {input.SNPS} \
            --allow-extra-chr \
            --recode \
            --out {params.output_prefix}) > {log} 2>&1"