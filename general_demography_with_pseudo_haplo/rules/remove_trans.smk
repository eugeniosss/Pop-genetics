rule remove_trans:
    input:
        FRQ="plink/"+config["subsets"]+".frq",
    output:
        FRQ=temp("plink/"+config["subsets"]+".snps"),
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
        BED="plink/"+config["subsets"]+".bed",
        BIM="plink/"+config["subsets"]+".bim",
        FAM="plink/"+config["subsets"]+".fam",
        SNPS="plink/"+config["subsets"]+".snps",
    output:
        MAP=temp("plink/"+config["subsets"]+"_no_trans.map"),
        NOSEX=temp("plink/"+config["subsets"]+"_no_trans.nosex"),
        PED=temp("plink/"+config["subsets"]+"_no_trans.ped"),
        LOG="plink/"+config["subsets"]+"_no_trans.log"
    log:
        "plink/"+config["subsets"]+"_no_trans_snake.log"
    conda:
        config["dir"] + "envs/plink.yml"
    params:
        input_prefix="plink/"+config["subsets"],
        output_prefix="plink/"+config["subsets"]+"_no_trans",
    shell:
        "(plink \
            --bfile {params.input_prefix} \
            --extract {input.SNPS} \
            --allow-extra-chr \
            --recode \
            --out {params.output_prefix}) > {log} 2>&1"