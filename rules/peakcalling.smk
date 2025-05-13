# rules/peakcalling.smk
# Peak calling rules for ChIP-seq data

rule macs3_callpeak_raw:
    input:
        "bam/{sample}_filtered.bam"
    output:
        raw_summits=temp("peaks/raw/{sample}.bed")
    params:
        prefix=lambda wc: wc.sample,
        genome_size=MACS3_GENOME,
        qvalue=MACS3_QVALUE
    shell:
        """
        macs3 callpeak -t {input} -f BAM -g {params.genome_size} \
            -n {params.prefix} -q {params.qvalue} --outdir peaks
        mv peaks/{params.prefix}_summits.bed {output.raw_summits}
        """

rule filter_summits_chromosomes:
    input:
        raw_summits="peaks/raw/{sample}.bed"
    output:
        summits_bed="peaks/{sample}_summits.bed"
    shell:
        """
        awk '$1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/' {input.raw_summits} > {output.summits_bed}
        """