


rule modify_BIM_all_to_1:
    input:
        BIM="plink/"+config["subsets"]+"_no_trans_pseudo_LD.bim",
    output:
        BIM="proje_pca/"+config["subsets"]+"_all_1.bim",
    shell:
        """
        awk '{{ $1="1"; }} 1' < {input.BIM} > {output.BIM}
        """

rule prepare_parfile_convertf:
    input:
        BED="plink/"+config["subsets"]+"_no_trans_pseudo_LD.bed",
        BIM="proje_pca/"+config["subsets"]+"_all_1.bim",
        FAM="plink/"+config["subsets"]+"_no_trans_pseudo_LD.fam",
    output:
        "proje_pca/"+config["subsets"]+"_bed_to_eigen.par",
    params:
        OUT_PREFIX="proje_pca/"+config["subsets"]
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
        "proje_pca/"+config["subsets"]+"_bed_to_eigen.par",
    output:
        EIGENSTRATGENO="proje_pca/"+config["subsets"]+".eigenstratgeno",
        SNP="proje_pca/"+config["subsets"]+".snp",
        IND="proje_pca/"+config["subsets"]+".ind",
    conda:
        config["dir"] + "envs/eigensoft.yml"
    log:
        "proje_pca/"+config["subsets"]+"_convertf.log"
    shell:
        "(convertf -p {input}) > {log} 2>&1"

rule add_categories_to_ind_file:
    input:
        IND="proje_pca/"+config["subsets"]+".ind",
        CORRESPONDENCIES=config["subsets"]+".txt",
    output:
        "proje_pca/"+config["subsets"]+"_updated.ind",
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
        EIGENSTRATGENO="proje_pca/"+config["subsets"]+".eigenstratgeno",
        SNP="proje_pca/"+config["subsets"]+".snp",
        IND="proje_pca/"+config["subsets"]+"_updated.ind",
    output:
        "proje_pca/"+config["subsets"]+"_smartpca.par",
    params:
        OUT_PREFIX="proje_pca/"+config["subsets"],
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
        "proje_pca/"+config["subsets"]+"_smartpca.par",
    output:
        EVEC="proje_pca/"+config["subsets"]+".evec",
        EVAL="proje_pca/"+config["subsets"]+".eval",
    conda:
        config["dir"] + "envs/eigensoft.yml"
    log:
        "proje_pca/"+config["subsets"]+"_smartpca.log"
    shell:
        "(smartpca -p {input}) > {log} 2>&1"
