rule randomize_haplotypes:
    input:
        "plink_{subset}/{subset}_no_trans.ped",
    output:
        temp("plink_{subset}/{subset}_tmp.txt"),
    run:
        import csv
        import random
        import itertools as it

        in_file = input[0]
        out_file = output[0]

        with open(out_file, 'w') as o_f:
            with open(in_file, 'r') as in_f:
                for line in csv.reader(in_f, dialect="excel", delimiter=' '):
                    out_head = line[0:6]
                    in_seq = line[6:]
                    in_geno = list(zip(*[iter(in_seq)]*2))
                    out_hap = [random.choice(el) for el in in_geno]
                    out_geno = list(zip(out_hap, out_hap))
                    out_seq = list(it.chain.from_iterable(out_geno))
                    out_line = " ".join(out_head + out_seq)
                    o_f.write(out_line + "\n")

rule overwrite_ped_pseudo_haploid:
    input:
        PED="plink_{subset}/{subset}_tmp.txt",
        MAP="plink_{subset}/{subset}_no_trans.map",
        NOSEX="plink_{subset}/{subset}_no_trans.nosex",
        LOG="plink_{subset}/{subset}_no_trans.log"
    output:
        PED=temp("plink_{subset}/{subset}_no_trans_mod.ped"),
        MAP=temp("plink_{subset}/{subset}_no_trans_mod.map"),
        NOSEX=temp("plink_{subset}/{subset}_no_trans_mod.nosex"),
        LOG="plink_{subset}/{subset}_no_trans_mod.log"
    shell:
        """
        mv {input.PED} {output.PED} 
        mv {input.LOG} {output.LOG} 
        mv {input.MAP} {output.MAP} 
        mv {input.NOSEX} {output.NOSEX} 
        """

rule pseudo_haploid_with_plink:
    input:
        PED="plink_{subset}/{subset}_no_trans_mod.ped",
        MAP="plink_{subset}/{subset}_no_trans_mod.map",
        NOSEX="plink_{subset}/{subset}_no_trans_mod.nosex",
        LOG="plink_{subset}/{subset}_no_trans_mod.log"
    output:
        BED="plink_{subset}/{subset}_no_trans_pseudo.bed",
        BIM="plink_{subset}/{subset}_no_trans_pseudo.bim",
        FAM="plink_{subset}/{subset}_no_trans_pseudo.fam",
        NOSEX="plink_{subset}/{subset}_no_trans_pseudo.nosex",
    params:
        input_prefix="plink_{subset}/{subset}_no_trans_mod",
        output_prefix="plink_{subset}/{subset}_no_trans_pseudo",
    conda:
        config["dir"] + "envs/plink.yml"
    log:
        "plink_{subset}/{subset}_no_trans_pseudo_snake.log",
    shell:
        "(plink \
        --file {params.input_prefix} \
        --allow-extra-chr \
        --maf 0.001 \
        --make-bed \
        --out {params.output_prefix}) > {log} 2>&1"