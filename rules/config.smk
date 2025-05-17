# rules/config.smk
# Configuration file for the Snakemake pipeline

configfile: "config.yaml"

import pandas as pd

# Read TSV, skip comment lines, and fill missing inputs with "NA"
samples_df = pd.read_csv(config["samples_tsv"], sep="\t", comment="#", skip_blank_lines=True).fillna("NA")

SAMPLES = samples_df["Sample"].tolist()
INPUTS = samples_df["Input"].tolist()
SAMPLE_TO_INPUT = dict(zip(SAMPLES, INPUTS))

ALL_GSMS = sorted({
    s.strip()                        # trim whitespace
    for s in SAMPLES + INPUTS
    if s and s.strip().upper() != "NA"
})


# Load config variables
BOWTIE2_INDEX = config["bowtie2_index"]
CHROM_SIZES = config["chrom_sizes"]
THREADS = config["threads"]
MACS3_GENOME = config["macs3_genome_size"]
MACS3_QVALUE = config["macs3_qvalue"]
FASTA_FILE = config["fasta_file"]
BPNET_CONFIG_GIN = config["bpnet_config"]
BPNET_CONTAINER = config["apptainer_sif"]
COMET_PROJECT = config["comet_project"]