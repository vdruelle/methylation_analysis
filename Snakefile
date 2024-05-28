# Pipeline for the analysis of the .bam generated using modified basecalling
import pathlib
import numpy as np
import os


#  check that config file is provided
def check_configfile():
    try:
        msg = f"""
----------- loading config file -----------
--- run config:
{config["run_config"]}
-------------------------------------------
"""
        print(msg)
    except:
        raise Exception(
            "config file not specified. Please specify with --configfile flag."
        )


def extract_record_names(bam_file):
    """Given a fasta file returns the list of record names.
    Each name is split at the first space occurrence.
    Checks that there are no duplicates."""
    records = []
    with open(bam_file, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                rec_name = line[1:].strip()
                assert (
                    " " not in rec_name
                ), f"space character not allowed in record name: {rec_name} from {bam_file}"
                records.append(rec_name)
    assert len(records) == len(
        np.unique(records)
    ), f"warning, duplicated record name in {bam_file}"
    return records


def parse_pileup_config(rc):
    """parser for pileup part of the config file."""
    # extract fasta record names for each reference
    records = {}
    for ref_name in rc["pileups"]:
        ref_file = rc["input"] + "/references/" + ref_name + ".fa"
        records[ref_name] = extract_record_names(ref_file)
    rc["ref_records"] = records

    return rc


def check_input_files_exist(run_config):
    input_files = []

    # Collect read files
    for sample_id in run_config["pileups"].values():
        for sid in sample_id:
            read_file = f"{run_config['input']}/reads/{sid}.bam"
            input_files.append(read_file)

    # Collect reference files
    for ref_id in run_config["pileups"].keys():
        ref_file = f"{run_config['input']}/references/{ref_id}.fa"
        input_files.append(ref_file)

    # Collect region files
    for ref_id in run_config["pileups"].keys():
        region_file = f"{run_config['input']}/regions/{ref_id}.bed"
        input_files.append(region_file)

    # Check if files exist
    missing_files = [f for f in input_files if not os.path.exists(f)]
    if missing_files:
        raise FileNotFoundError(
            f"The following input files are missing: {', '.join(missing_files)}"
        )
    else:
        print("All input files are present.")


check_configfile()
run_config = parse_pileup_config(config["run_config"])
check_input_files_exist(run_config)
pathlib.Path("log").mkdir(exist_ok=True)

INPUT_DIR = run_config["input"].removesuffix("/")
OUTPUT_DIR = run_config["output"].removesuffix("/")
READS_DIR = INPUT_DIR + "/reads"
REFERENCES_DIR = INPUT_DIR + "/references"
REGION_DIR = INPUT_DIR + "/regions"
ALIGNMENTS_DIR = OUTPUT_DIR + "/alignments"
PILEUP_DIR = OUTPUT_DIR + "/pileup"
FOCUS_REGIONS_DIR = OUTPUT_DIR + "/focus_regions"

DORADO_BIN = "softwares/dorado-0.7.0-linux-x64/bin/dorado"
MODKIT_BIN = "softwares/modkit/modkit"


rule all:
    input:
        [
            expand(
                FOCUS_REGIONS_DIR + "/{ref_id}/{sample_id}.bed",
                ref_id=ref_id,
                sample_id=sample_id,
            )
            for ref_id, sample_id in run_config["pileups"].items()
        ],


rule align:
    message:
        "Aligning reads {input.reads} to the reference {input.ref}."
    input:
        reads=READS_DIR + "/{sample_id}.bam",
        ref=REFERENCES_DIR + "/{ref_id}.fa",
    output:
        bam=ALIGNMENTS_DIR + "/{ref_id}/{sample_id}.bam",
        bam_sorted=ALIGNMENTS_DIR + "/{ref_id}/{sample_id}_sorted.bam",
        bai=ALIGNMENTS_DIR + "/{ref_id}/{sample_id}_sorted.bam.bai",
    params:
        dorado=DORADO_BIN,
    threads: 4
    conda:
        "methylation_analysis"
    shell:
        """
        ./{params.dorado} aligner {input.ref} {input.reads} -t {threads} > {output.bam}
        samtools sort {output.bam} > {output.bam_sorted}
        samtools index {output.bam_sorted}
        """


rule pileup:
    message:
        "Performing pileup of alignment {input.bam} using modkit."
    input:
        bam=rules.align.output.bam_sorted,
        ref=rules.align.input.ref,
    output:
        pileup=PILEUP_DIR + "/{ref_id}/{sample_id}.pileup",
    params:
        modkit=MODKIT_BIN,
    threads: 4
    shell:
        """
        ./{params.modkit} pileup {input.bam} {output.pileup} -r {input.ref} -t {threads}
        """


rule intersect:
    message:
        "Extracting regions defined in {input.regions} from {input.pileup} using bedtools intersect."
    input:
        pileup=rules.pileup.output.pileup,
        regions=REGION_DIR + "/{ref_id}.bed",
    output:
        intersect=FOCUS_REGIONS_DIR + "/{ref_id}/{sample_id}.bed",
    conda:
        "methylation_analysis"
    shell:
        """
        bedtools intersect -a {input.pileup} -b {input.regions} > {output.intersect}
        """
