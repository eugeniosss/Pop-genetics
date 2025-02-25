
rule sub_sample_vcf:
	input:
		vcf="{software}_{dataset}/variants.vcf.gz",
		tbi="{software}_{dataset}/variants.vcf.gz.tbi"
	output:
		vcf="{software}_{dataset}/variants_sub.vcf.gz",
		tbi="{software}_{dataset}/variants_sub.vcf.gz.tbi"
	conda:
		config["dir"] + "envs/geno_callers_env.yml"
	threads: 1
	params:
		percentage=config["options"]["sub_sample_vcf_for_stats"]
	log:
		"{software}_{dataset}/variants_sub.log"
	shell:
		"""
    	(bcftools view {input.vcf} | vcfrandomsample -r {params.percentage} | bgzip -c > {output.vcf} &&
    	tabix -p vcf {output.vcf}) > {log} 2>&1
    	"""

rule calculate_stats_vcf:
	input:
		"{software}_{dataset}/variants_sub.vcf.gz",
	output:
		lqual="{software}_{dataset}/stats/site_quality.lqual",
		ldepth_mean="{software}_{dataset}/stats/mean_depth_per_site.ldepth.mean",
		lmiss="{software}_{dataset}/stats/missing_data_per_site.lmiss",
		frq="{software}_{dataset}/stats/allele_freq2.frq",
		idepth="{software}_{dataset}/stats/mean_depth_per_individual.idepth",
		imiss="{software}_{dataset}/stats/missing_data_per_indv.imiss"
	threads: 1
	params:
		output_folder="{software}_{dataset}/stats"
	log:
		"{software}_{dataset}/stats/calculating.log"
	shell:
		"""
   		# Calculate allele frequency
		(vcftools --gzvcf {input} --freq2 --out {params.output_folder}/allele_freq2 --min-alleles 2 --max-alleles 2

		# Mean depth per individual
		vcftools --gzvcf {input} --depth --out {params.output_folder}/mean_depth_per_individual
		
		# Mean depth per site
		vcftools --gzvcf {input} --site-mean-depth --out {params.output_folder}/mean_depth_per_site
		
		# Calculating site quality
		vcftools --gzvcf {input} --site-quality --out {params.output_folder}/site_quality
		
		# Proportion of Missing Data per indv
		vcftools --gzvcf {input} --missing-indv --out {params.output_folder}/missing_data_per_indv
		
		# Missing data per site
		vcftools --gzvcf {input} --missing-site --out {params.output_folder}/missing_data_per_site) > {log} 2>&1

    	"""

rule create_ind_plot:
	input:
		lqual="{software}_{dataset}/stats/site_quality.lqual",
		ldepth_mean="{software}_{dataset}/stats/mean_depth_per_site.ldepth.mean",
		lmiss="{software}_{dataset}/stats/missing_data_per_site.lmiss",
		frq="{software}_{dataset}/stats/allele_freq2.frq",
		idepth="{software}_{dataset}/stats/mean_depth_per_individual.idepth",
		imiss="{software}_{dataset}/stats/missing_data_per_indv.imiss"
	output:
		"{software}_{dataset}/stats/{software}_{dataset}_stats_summarized.pdf"
	conda:
		config["dir"] + "envs/R_general.yaml"
	params:
		input_dir="{software}_{dataset}/stats/",
		script=config["dir"]+"rules/summarize_stats_vcf.R"
	log:
		"{software}_{dataset}/stats/R_plotting.log"
	shell:
		"""
		(Rscript {params.script} {params.input_dir} {output}) > {log} 2>&1
    	"""

rule create_bed_files_for_positions_matrix:
	input:
		vcf="{software}_{dataset}/variants.vcf.gz",
		tbi="{software}_{dataset}/variants.vcf.gz.tbi"
	output:
		"{software}_{dataset}/variants.bed"
	conda:
		config["dir"] + "envs/geno_callers_env.yml"
	log:
		"{software}_{dataset}/creating_bed_file_soft.log"
	shell:
		"""
		echo {wildcards.software}_{wildcards.dataset} > {output}
		bcftools view -H {input.vcf} | awk '{{print $1":"$2"-"$2+1}}' >> {output}
		"""

rule create_positions_matrix_per_dataset:
	input:
		beds=expand("{software}_{{dataset}}/variants.bed", software=config['softwares'].split(" ")),
	output:
		"{dataset}_positions_matrix_soft.txt"
	log:
		"{dataset}_creating_intersect_matrix.log"
	shell:
		"""
		#!/usr/bin/bash
	
		paste {input.beds} > {output} 
		"""

rule create_multi_plot_per_dataset:
	input:
		lqual=expand("{software}_{{dataset}}/stats/site_quality.lqual", software=config['softwares'].split(" ")),
		ldepth_mean=expand("{software}_{{dataset}}/stats/mean_depth_per_site.ldepth.mean", software=config['softwares'].split(" ")),
		lmiss=expand("{software}_{{dataset}}/stats/missing_data_per_site.lmiss", software=config['softwares'].split(" ")),
		frq=expand("{software}_{{dataset}}/stats/allele_freq2.frq", software=config['softwares'].split(" ")),
		idepth=expand("{software}_{{dataset}}/stats/mean_depth_per_individual.idepth", software=config['softwares'].split(" ")),
		imiss=expand("{software}_{{dataset}}/stats/missing_data_per_indv.imiss", software=config['softwares'].split(" ")),
		intersects="{dataset}_positions_matrix_soft.txt"
	output:
		"{dataset}_stats_summarized_different_softwares.pdf"
	conda:
		config["dir"] + "envs/R_general.yaml"
	params:
		input_dir=",".join(expand("{software}_{{dataset}}/stats/", software=config['softwares'].split(" "))),
		script=config["dir"]+"rules/summarize_multiple_stats_vcf_per_dataset.R"
	log:
		"{dataset}_R_plotting.log"
	shell:
		"""
		(Rscript {params.script} -f {params.input_dir} -m {input.intersects} -o {output}) > {log} 2>&1
    	"""

rule create_positions_matrix_per_software:
	input:
		beds=expand("{{software}}_{dataset}/variants.bed", dataset=config['datasets'].split(" ")),
	output:
		"{software}_positions_matrix_data.txt"
	log:
		"{software}_creating_intersect_matrix.log"
	shell:
		"""
		#!/usr/bin/bash
	
		paste {input.beds} > {output} 
		"""

rule create_multi_plot_per_software:
	input:
		lqual=expand("{{software}}_{dataset}/stats/site_quality.lqual", dataset=config['datasets'].split(" ")),
		ldepth_mean=expand("{{software}}_{dataset}/stats/mean_depth_per_site.ldepth.mean", dataset=config['datasets'].split(" ")),
		lmiss=expand("{{software}}_{dataset}/stats/missing_data_per_site.lmiss", dataset=config['datasets'].split(" ")),
		frq=expand("{{software}}_{dataset}/stats/allele_freq2.frq", dataset=config['datasets'].split(" ")),
		idepth=expand("{{software}}_{dataset}/stats/mean_depth_per_individual.idepth", dataset=config['datasets'].split(" ")),
		imiss=expand("{{software}}_{dataset}/stats/missing_data_per_indv.imiss", dataset=config['datasets'].split(" ")),
		intersects="{software}_positions_matrix_data.txt"
	output:
		"{software}_stats_summarized_different_datasets.pdf"
	conda:
		config["dir"] + "envs/R_general.yaml"
	params:
		input_dir=",".join(expand("{{software}}_{dataset}/stats/", dataset=config['datasets'].split(" "))),
		script=config["dir"]+"rules/summarize_multiple_stats_vcf_per_software.R"
	log:
		"{software}_R_plotting.log"
	shell:
		"""
		(Rscript {params.script} -f {params.input_dir} -m {input.intersects} -o {output}) > {log} 2>&1
    	"""