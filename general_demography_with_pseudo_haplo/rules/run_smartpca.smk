


rule modify_BIM_all_to_1:
    input:
        BIM="plink_{subset}/{subset}_no_trans_pseudo_LD.bim",
    output:
        BIM="ProjePca_{subset}/{subset}_all_1.bim",
    shell:
        """
        awk '{{ $1="1"; }} 1' < {input.BIM} > {output.BIM}
        """

rule prepare_parfile_convertf:
    input:
        BED="plink_{subset}/{subset}_no_trans_pseudo_LD.bed",
        BIM="ProjePca_{subset}/{subset}_all_1.bim",
        FAM="plink_{subset}/{subset}_no_trans_pseudo_LD.fam",
    output:
        "ProjePca_{subset}/{subset}_bed_to_eigen.par",
    params:
        OUT_PREFIX="ProjePca_{subset}/{subset}"
    shell:
        """
        echo genotypename: {input.BED} > {output}
        echo snpname: {input.BIM} >> {output}
        echo indivname: {input.FAM} >> {output}
        echo outputformat: EIGENSTRAT >> {output}
        echo genotypeoutname: {params.OUT_PREFIX}.eigenstratgeno >> {output}
        echo snpoutname: {params.OUT_PREFIX}.snp >> {output}
        echo indivoutname: {params.OUT_PREFIX}.ind >> {output}
        echo familynames: NO >> {output}
        echo pordercheck: NO >> {output}
        """

rule run_convertf:
    input:
        "ProjePca_{subset}/{subset}_bed_to_eigen.par",
    output:
        EIGENSTRATGENO="ProjePca_{subset}/{subset}.eigenstratgeno",
        SNP="ProjePca_{subset}/{subset}.snp",
        IND="ProjePca_{subset}/{subset}.ind",
    conda:
        config["dir"] + "envs/eigensoft.yml"
    log:
        "ProjePca_{subset}/{subset}_convertf.log"
    shell:
        "(convertf -p {input}) > {log} 2>&1"

rule add_categories_to_ind_file:
    input:
        IND="ProjePca_{subset}/{subset}.ind",
        CORRESPONDENCIES="{subset}.txt",
    output:
        "ProjePca_{subset}/{subset}_updated.ind",
    shell:
        """       
        cp {input.IND} {output}

        while IFS= read -r line || [[ -n "$line" ]]; do

            sample=$(echo $line | cut -f 1 -d " ")

            group=$(echo $line | cut -f 2 -d " ")

            sed -i "/^[ ]*${{sample}}/s/???/${{group}}/g" {output}

        done < {input.CORRESPONDENCIES}
        """

rule prepare_parfile_smartpca:
    input:
        EIGENSTRATGENO="ProjePca_{subset}/{subset}.eigenstratgeno",
        SNP="ProjePca_{subset}/{subset}.snp",
        IND="ProjePca_{subset}/{subset}_updated.ind",
    output:
        "ProjePca_{subset}/{subset}_smartpca.par",
    params:
        OUT_PREFIX="ProjePca_{subset}/{subset}",
    shell:
        """
            echo genotypename:      {input.EIGENSTRATGENO} > {output}
            echo snpname:           {input.SNP} >> {output}
            echo indivname:         {input.IND} >> {output}
            echo evecoutname:       {params.OUT_PREFIX}.evec >> {output}
            echo evaloutname:       {params.OUT_PREFIX}.eval >> {output}
            echo poplistname:       times.txt >> {output}
            echo outliermode: 2 >> {output}
            echo lsqproject:         YES >> {output}

            echo Modern > times.txt
        """

rule run_smartpca:
    input:
        "ProjePca_{subset}/{subset}_smartpca.par",
    output:
        EVEC="ProjePca_{subset}/{subset}.evec",
        EVAL="ProjePca_{subset}/{subset}.eval",
    conda:
        config["dir"] + "envs/eigensoft.yml"
    log:
        "ProjePca_{subset}/{subset}_smartpca.log"
    shell:
        "(smartpca -p {input}) > {log} 2>&1"
