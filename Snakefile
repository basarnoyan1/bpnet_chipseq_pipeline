include: "rules/config.smk"
include: "rules/indexing.smk"
include: "rules/preprocessing.smk"
include: "rules/peakcalling.smk"
include: "rules/bpnet.smk"

rule all:
    input:
        expand("sra_id/{sample}.sra_id", sample=ALL_GSMS),
        expand("bam/{sample}_filtered.bam", sample=ALL_GSMS),
        expand("bam/{sample}_filtered.bam.bai", sample=ALL_GSMS),
        expand("peaks/{sample}_summits.bed", sample=SAMPLES),
        expand("bigwig/{sample}_plus.bw", sample=ALL_GSMS),
        expand("bigwig/{sample}_minus.bw", sample=ALL_GSMS),
        expand("dataspec/{sample}_dataspec.yml", sample=SAMPLES),
        expand("bpnet_training/{sample}_complete", sample=SAMPLES),
        expand("contrib/{sample}.h5", sample=SAMPLES),
        expand("contrib/{sample}_null.h5", sample=SAMPLES),
        expand("modisco/{sample}/modisco_done", sample=SAMPLES),