##################################
# rules/indexing.smk
# Build Bowtie2 index files from reference FASTA
##################################

# Assumes FASTA_FILE and BOWTIE2_INDEX are defined in config
rule build_bowtie2_index:
    """
    Build Bowtie2 index from reference FASTA.
    """
    input:
        fasta=FASTA_FILE
    output:
        # pick one of the six .bt2 shards as the marker
        f"{BOWTIE2_INDEX}.1.bt2"
    params:
        prefix=BOWTIE2_INDEX
    threads: 4
    log:
        "logs/indexing/bowtie2-build.log"
    shell:
        """
        mkdir -p $(dirname {output})
        bowtie2-build --threads {threads} {input.fasta} {params.prefix} \
            &> {log}
        """
