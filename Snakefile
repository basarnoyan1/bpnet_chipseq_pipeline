include: "rules/config.smk"
include: "rules/preprocessing.smk"
include: "rules/peakcalling.smk"
include: "rules/bpnet.smk"

rule all:
    input:
        expand("sra_id/{sample}.sra_id", sample=ALL_SRAS),
        expand("bam/{sample}_filtered.bam", sample=ALL_SRAS),
        expand("bam/{sample}_filtered.bam.bai", sample=ALL_SRAS),
        expand("peaks/{sample}_summits.bed", sample=SAMPLES),
        expand("bigwig/{sample}_plus.bw", sample=ALL_SRAS),
        expand("bigwig/{sample}_minus.bw", sample=ALL_SRAS),
        expand("dataspec/{sample}_dataspec.yml", sample=SAMPLES),
        expand("bpnet_training/{sample}_complete", sample=SAMPLES),
        expand("contrib/{sample}.h5", sample=SAMPLES),
        expand("contrib/{sample}_null.h5", sample=SAMPLES),
        expand("modisco/{sample}/chip_nexus_done", sample=SAMPLES)