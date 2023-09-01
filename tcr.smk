rule all:
	input:
		expand("clean_data/{sample}_clean.1.fq.gz", sample=['A10','A11','A1','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A2','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A3','A32','A33','A34','A35','A36','A37','A38','A39','A40','A4','A5','A6','A7','A8','A9']),
		expand("MiXCR/{sample}_.clonotypes.ALL.txt", sample=['A10','A11','A1','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A2','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A3','A32','A33','A34','A35','A36','A37','A38','A39','A40','A4','A5','A6','A7','A8','A9']),
		expand("clonotype/{sample}.TRB.tsv", sample=['A10','A11','A1','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A2','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A3','A32','A33','A34','A35','A36','A37','A38','A39','A40','A4','A5','A6','A7','A8','A9'])


rule fastp:
	input:
		fq1="raw_data/{sample}_1.fq.gz",
		fq2="raw_data/{sample}_2.fq.gz",
	output:
		fq1="clean_data/{sample}_clean.1.fq.gz",
		fq2="clean_data/{sample}_clean.2.fq.gz",
	threads:
		workflow.cores * 0.5
	log:
		"logs/{sample}.fastp.log"
	shell:
		"""
		/home/haoq/biosoftware/fastp -i {input.fq1} -o {output.fq1} \
			-I {input.fq2} -O {output.fq2} \
			-3 \
			-q 20 -u 30 \
			--length_required=25 \
			--thread {threads} \
			--compression=6 \
			-j clean_data/{wildcards.sample}.json -h clean_data/{wildcards.sample}.html 2> {log}
		"""



rule Mixcr:
	input:
		fq1="clean_data/{sample}_clean.1.fq.gz",
		fq2="clean_data/{sample}_clean.2.fq.gz"
	output:
		"MiXCR/{sample}_.clonotypes.ALL.txt"
	log:
		"logs/{sample}_.Mixcr.log"
	shell:
		"""
		mixcr analyze amplicon \
            --species hs \
            --starting-material rna \
            --5-end v-primers \
            --3-end j-primers \
            --adapters no-adapters \
            {input.fq1} {input.fq2} MiXCR/{wildcards.sample}_  2> {log}
		"""

rule clonotype:
    input:
        "MiXCR/{sample}_.clns"
    output:
        "clonotype/{sample}.TRB.tsv"
    shell:
        '''
        mixcr exportClones --chains TRB {input} \
			-o -t \
            clonotype/{wildcards.sample}.TRB.tsv
        '''

