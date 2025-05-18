# rules/preprocessing.smk
# Requires globals: THREADS, BOWTIE2_INDEX, CHROM_SIZES
############################################################

##############################
# 1. fetch_sra  (pysradb lookup + download)
##############################
rule fetch_sra:
    """
    Resolve SRR run IDs for each GSM using pysradb, then download the .sra files.
    """
    output:
        ids    = "sra_id/{sample}.sra_id",
        marker = touch("sra/{sample}/download.done")
    params:
        gsm = lambda wc: wc.sample
    threads: 1
    run:
        import subprocess, pathlib, sys

        gsm   = params.gsm
        ids_f = output.ids
        done  = output.marker

        # Skip placeholder NA
        if gsm == "NA":
            pathlib.Path(ids_f).write_text("NA\n")
            pathlib.Path(done).touch()
            return

        # Lookup with pysradb
        try:
            out = subprocess.check_output(
                ["pysradb", "gsm-to-srr", gsm],
                text=True
            ).strip().splitlines()
        except FileNotFoundError:
            sys.exit("ERROR: pysradb not found (conda install -c bioconda pysradb).")
        except subprocess.CalledProcessError as e:
            sys.exit(f"ERROR: pysradb failed for {gsm}: {e}")

        if len(out) <= 1:
            sys.exit(f"ERROR: No SRR runs found for {gsm}")

        # Parse SRR column (skip header)
        srrs = [line.split()[1] for line in out[1:]]
        pathlib.Path("sra_id").mkdir(exist_ok=True)
        with open(ids_f, "w") as fh:
            for srr in srrs:
                fh.write(srr + "\n")

        # Download each SRR
        dest = pathlib.Path(f"sra/{gsm}")
        dest.mkdir(parents=True, exist_ok=True)
        for srr in srrs:
            subprocess.check_call(
                ["prefetch", srr, "--output-file", str(dest / f"{srr}.sra")]
            )

        pathlib.Path(done).touch()


##############################
# 2. sra → per-run FASTQs
##############################
rule sra2fastq:
    input:
        marker="sra/{sample}/download.done",
        ids   ="sra_id/{sample}.sra_id"
    output:
        touch("fastq/{sample}/fastq.done")
    threads: THREADS
    shell:
        r"""
        if grep -q "^NA$" {input.ids}; then
            mkdir -p fastq/{wildcards.sample}; touch {output}; exit 0
        fi
        mkdir -p fastq/{wildcards.sample}
        for srr in $(cat {input.ids}); do
            if [ ! -f fastq/{wildcards.sample}/${{srr}}_1.fastq.gz ]; then
                python parallel-fastq-dump.py -s sra/{wildcards.sample}/${{srr}}.sra \
                   --outdir fastq/{wildcards.sample} --split-files --gzip --threads {threads}
            fi
        done
        touch {output}
        """

rule align_and_merge:
    """
    Align all FASTQs for a sample (in fastq/{sample}/) and merge into bam/{sample}_merged.bam.
    """
    input:
        fastq_done="fastq/{sample}/fastq.done"
    output:
        merged_bam=temp("bam/{sample}_merged.bam")
    params:
        index=BOWTIE2_INDEX
    threads: THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p bam/runs/{wildcards.sample}

        for fq1 in fastq/{wildcards.sample}/*_1.fastq.gz; do
            run=$(basename "$fq1" _1.fastq.gz)
            fq2="fastq/{wildcards.sample}/${{run}}_2.fastq.gz"
            sam="bam/runs/{wildcards.sample}/${{run}}.sam"
            sorted="bam/runs/{wildcards.sample}/${{run}}_sorted.bam"

            # 1) align to SAM
            if [ -f "$fq2" ]; then
                bowtie2 -x {params.index} -1 "$fq1" -2 "$fq2" \
                    -p {threads} --quiet -S "$sam"
            else
                bowtie2 -x {params.index} -U "$fq1" \
                    -p {threads} --quiet -S "$sam"
            fi

            # 2) convert & sort
            samtools view -bS -q 30 -@ {threads} "$sam" \
                | samtools sort -@ {threads} -o "$sorted"

            rm "$sam"
        done

        # 3) merge all per-run BAMs
        samtools merge -@ {threads} "{output.merged_bam}" \
            bam/runs/{wildcards.sample}/*_sorted.bam
        """


rule filter_chromosomes:
    input:
        merged="bam/{sample}_merged.bam"
    output:
        "bam/{sample}_filtered.bam"
    threads: THREADS
    shell:
        """
        samtools view -@ {threads} -h {input.merged} \
         | awk '/^@/ || $3 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/' \
         | samtools view -b -@ {threads} -o {output}
        """


rule index_bam:
    input:  "bam/{sample}_filtered.bam"
    output: "bam/{sample}_filtered.bam.bai"
    shell:  "samtools index {input}"


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
        r"""
        bedtools genomecov -5 -bg -strand + -ibam {input} \
          | sort -k1,1 -k2,2n > tmp_{wildcards.sample}_plus.bedGraph
        bedGraphToBigWig \
          tmp_{wildcards.sample}_plus.bedGraph {params.sizes} {output.plus_bw}

        bedtools genomecov -5 -bg -strand - -ibam {input} \
          | sort -k1,1 -k2,2n > tmp_{wildcards.sample}_minus.bedGraph
        bedGraphToBigWig \
          tmp_{wildcards.sample}_minus.bedGraph {params.sizes} {output.minus_bw}

        rm tmp_{wildcards.sample}_*.bedGraph
        """

