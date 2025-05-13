# rules/preprocessing.smk
# This file contains rules for preprocessing sequencing data, including downloading SRA files,
# converting SRA to FASTQ, aligning reads with Bowtie2, filtering chromosomes, and generating bigWig files.

rule retrieve_sra_code:
    output:
        "sra_id/{sample}.sra_id"
    params:
        query=lambda wc: wc.sample
    threads: 1
    retries: 5
    shell:
        """
        echo "[INFO] Retrieving SRA ID for sample: {params.query}"
        esearch -db sra -query "{params.query}" | \
        efetch -format runinfo | grep -v "^Run" | cut -d ',' -f1 > {output}
        if [ ! -s {output} ]; then
            echo "[ERROR] No SRA Run ID found for {params.query}" >&2
            exit 1
        fi
        """

rule download_sra:
    input:
        "sra_id/{sample}.sra_id"
    output:
        "sra/{sample}.sra"
    threads: THREADS
    retries: 5
    shell:
        "prefetch $(cat {input}) --output-file {output}"

rule sra2fastq:
    input:
        "sra/{sample}.sra"
    output:
        "fastq/{sample}_1.fastq.gz" 
    threads: THREADS
    shell:
        """
        python parallel-fastq-dump.py -s {input} \
            --outdir fastq \
            --split-files \
            --gzip \
            --threads {threads}
        """

rule bowtie2_align:
    input:
        fastq1="fastq/{sample}_1.fastq.gz",
        fastq2="fastq/{sample}_2.fastq.gz"
    output:
        bam="bam/tmp_{sample}_sorted.bam"
    params:
        index=BOWTIE2_INDEX
    threads: THREADS
    run:
        import os
        sam = f"tmp_{wildcards.sample}.sam"
        if os.path.exists(input.fastq2):
            bowtie_cmd = f"bowtie2 -x {params.index} -1 {input.fastq1} -2 {input.fastq2} -p {threads} --verbose -S {sam}"
        else:
            bowtie_cmd = f"bowtie2 -x {params.index} -U {input.fastq1} -p {threads} --verbose -S {sam}"
        shell(f"""
            {bowtie_cmd}
            samtools view -bS -q 30 -@ {threads} {sam} | samtools sort -@ {threads} -o {output.bam}
            rm {sam}
        """)

rule filter_chromosomes:
    input:
        "bam/tmp_{sample}_sorted.bam"
    output:
        "bam/{sample}_filtered.bam"
    threads: THREADS
    shell:
        """
        samtools view -@ {threads} -h {input} | \
        awk '/^@/ || $3 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/' | \
        samtools view -b -@ {threads} -o {output}
        """

rule index_bam:
    input:
        "bam/{sample}_filtered.bam"
    output:
        "bam/{sample}_filtered.bam.bai"
    threads: 1
    shell:
        "samtools index {input}"

rule bam2bigwig:
    input:
        "bam/{sample}_filtered.bam"
    output:
        plus_bw="bigwig/{sample}_plus.bw",
        minus_bw="bigwig/{sample}_minus.bw"
    params:
        sizes=CHROM_SIZES
    threads: THREADS
    shell:
        """
        bedtools genomecov -5 -bg -strand + -ibam {input} | sort -k1,1 -k2,2n > tmp_plus.bedGraph
        bedGraphToBigWig tmp_plus.bedGraph {params.sizes} {output.plus_bw}
        bedtools genomecov -5 -bg -strand - -ibam {input} | sort -k1,1 -k2,2n > tmp_minus.bedGraph
        bedGraphToBigWig tmp_minus.bedGraph {params.sizes} {output.minus_bw}
        rm tmp_*.bedGraph
        """