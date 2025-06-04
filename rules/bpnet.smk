# rules/bpnet.smk
# This file contains rules for training a BPNet model, generating contribution files, and running Modisco analysis.

rule generate_dataspec:
    output:
        "dataspec/{sample}_dataspec.yml"
    run:
        import yaml
        sample = wildcards.sample
        input_sample = SAMPLE_TO_INPUT[sample]
        dataspec = {
            'fasta_file': FASTA_FILE,
            'task_specs': {
                "AR": {
                    'tracks': [
                        f"bigwig/{sample}_plus.bw",
                        f"bigwig/{sample}_minus.bw"
                    ],
                    'peaks': f"peaks/{sample}_summits.bed"
                }
            }
        }
        if input_sample != "NA":
            dataspec['bias_specs'] = {
                "input": {
                    'tracks': [
                        f"bigwig/{input_sample}_plus.bw",
                        f"bigwig/{input_sample}_minus.bw"
                    ],
                    'tasks': ["AR"]
                }
            }
        with open(output[0], 'w') as f:
            yaml.dump(dataspec, f, default_flow_style=False)

rule bpnet_train:
    input:
        dataspec="dataspec/{sample}_dataspec.yml",
        config_gin=BPNET_CONFIG_GIN,
        plus_bw="bigwig/{sample}_plus.bw",
        minus_bw="bigwig/{sample}_minus.bw",
        peaks="peaks/{sample}_summits.bed"
    output:
        touch("bpnet_training/{sample}_complete")
    log: "logs/bpnet_training/{sample}.log"
    params:
        run_id=lambda wc: f"bpnet_{wc.sample}",
        container=BPNET_CONTAINER,
        comet_project=COMET_PROJECT
    threads: 1
    shell:
        """
        apptainer exec --nv {params.container} \
            env PATH="/opt/micromamba/envs/bpnet/bin:$PATH" bpnet train --config={input.config_gin} --run-id={params.run_id} \
            --memfrac-gpu=1 --cometml-project={params.comet_project} --overwrite \
            {input.dataspec} bpnet_training/{params.run_id} > {log}
        """

rule bpnet_contrib:
    input:
        "bpnet_training/{sample}_complete"
    output:
        contrib_file="contrib/{sample}.h5",
        contrib_null_file="contrib/{sample}_null.h5"
    log: "logs/bpnet_contrib/{sample}.log"
    params:
        container=BPNET_CONTAINER,
        model_dir=lambda wc: f"bpnet_training/bpnet_{wc.sample}/bpnet_{wc.sample}",
        memfrac_gpu=1
    shell:
        """
        apptainer exec --nv {params.container} \
            env PATH="/opt/micromamba/envs/bpnet/bin:$PATH" bpnet contrib {params.model_dir} --method=deeplift --memfrac-gpu={params.memfrac_gpu} {output.contrib_file} > {log}

        apptainer exec --nv {params.container} \
            env PATH="/opt/micromamba/envs/bpnet/bin:$PATH" bpnet contrib {params.model_dir} --method=deeplift --memfrac-gpu={params.memfrac_gpu} --shuffle-seq --max-regions 5000 {output.contrib_null_file} >> {log}
        """
rule bpnet_modisco_run:
    input:
        contrib="contrib/{sample}.h5",
        contrib_null="contrib/{sample}_null.h5"
    output:
        touch("modisco/{sample}/modisco_done")
    log:
        "logs/modisco/{sample}.log"
    params:
        container = BPNET_CONTAINER,
        modisco_dir = lambda wc: f"modisco/{wc.sample}/",
        gpu_frac    = 1
    threads: 1
    shell:
        r"""
        apptainer exec --nv {params.container} \
            env PATH="/opt/micromamba/envs/bpnet/bin:$PATH" bpnet modisco-run {input.contrib} \
              --null-contrib-file={input.contrib_null} \
              --contrib-wildcard=AR/counts/pre-act \
              --premade=modisco-50k \
              --only-task-regions {params.modisco_dir} \
              --overwrite \
              --memfrac-gpu={params.gpu_frac} \
            > {log}

        touch {output}
        """

