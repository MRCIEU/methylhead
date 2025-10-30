import os
import subprocess
import pandas


# ===================================================================
# File paths
# ===================================================================

# === normalize_path(path, base) =====================================
def normalize_path(path, base=None):
    if not base is None \
       and not os.path.isabs(path) \
       and os.path.expanduser(path) == path:
        base = normalize_path(base)
        path = f"{base}/{path}"
    else:     
        path = os.path.expanduser(path) # expand ~
        path = os.path.abspath(path)    # expand relative path
    path = os.path.realpath(path)       # resolve symbolic links
    return path

# === dirname(path): directory of file path =================
def dirname(path):
    return os.path.dirname(normalize_path(path))

# === ensure all config file paths are absolute =====================
for name,path in config['inputs'].items():
    config['inputs'][name] = normalize_path(path)

for name,path in config['paths'].items():
    config['paths'][name] = normalize_path(path)
    
# === create output directory =========================================
try:
    os.makedirs(config['paths']['output'])
except FileExistsError:
    pass
    
# === inputs:fasta--path to genome assembly fasta file ==============
config['inputs']['fasta'] = f"{config['paths']['genome']}/{config['parameters']['assembly']}.fa"

# === paths:scripts--directory of repository scripts ===================
config['paths']['scripts'] = workflow.basedir + "/scripts"


# =====================================================================
# Samples
# =====================================================================

# === all_samples: list of sample ids ===============================
samples = pandas.read_csv(config["inputs"]["samplesheet"])
all_samples = samples["sample_id"].unique().tolist()



# =====================================================================
# Apptainer 
# =====================================================================
    
# ==== pull container images ==========================================
def pull_image(url,image_file):
    if not os.path.exists(image_file):        
        subprocess.run(["apptainer","pull",image_file,url],check=True)
    return image_file

for container in config['containers'].values():
    pull_image(container['image_url'], container['image_file'])
    
# === apptainer_exec--command to run commands using apptainer ===
def apptainer_exec(analysis, cmd, extra_args="--writable-tmpfs"):
    bind_paths = [path for path in config['paths'].values()] \
        + [dirname(path) for path in config['inputs'].values()] \
        + [os.getcwd()]
    bind_paths = list(set(bind_paths))
    apptainer_args = (
        " ".join([ f"-B {path}:{path}" for path in bind_paths ]) +
        " --pwd " + os.getcwd())
    image = config['containers'][analysis]['image_file']
    return f"apptainer exec {apptainer_args} {extra_args} {image} /bin/bash -c \"{cmd}\""



# =====================================================================
# HPC resources
# =====================================================================

# === get_resources--hpc resources required =============================
def get_resources(resource_class="light", **overrides):
    resources = config["resources"][resource_class].copy()
    resources.update(overrides)
    return resources
        


# ===================================================================
# Rules
# ===================================================================

include: "rules/genome.smk"
include: "rules/merge_lanes.smk"
include: "rules/fastqc.smk"
include: "rules/trim.smk"
include: "rules/panel_intervals.smk"
include: "rules/bwa_align.smk"
include: "rules/sambamba.smk"
include: "rules/sort_bam.smk"
include: "rules/mark_duplicates.smk"
include: "rules/collect_hs_metrics.smk"
include: "rules/collect_mm_metrics.smk"
include: "rules/samtools_stats.smk"
include: "rules/multiqc.smk"
include: "rules/mbias.smk"
include: "rules/extract_cpg_bedgraph.smk"
include: "rules/meth_bedgraph.smk"
include: "rules/extract_cpg_methylation.smk"
include: "rules/methylation_matrix.smk"
include: "rules/illumina_matrix.smk"
include: "rules/estimate_cell_counts.smk"
include: "rules/dna_methylation_scores.smk"
include: "rules/bsmap_align.smk"
include: "rules/camda.smk"
include: "rules/camda_matrix.smk"
include: "rules/raw_read_counts.smk"
include: "rules/qc_report.smk"
include: "rules/test_associations.smk"


rule all:
    input:
        config['paths']['output'] + "/methylation_matrix/methylation-matrix.csv.gz",
        config['paths']['output'] + "/illumina_matrix/illumina-matrix.csv.gz",
        config['paths']['output'] + "/camda_matrix/camda-matrix.csv.gz",
        config['paths']['output'] + "/estimate_cell_counts/counts.csv",
        config['paths']['output'] + "/dna_methylation_scores/scores.csv",
        config['paths']['output'] + "/multiqc/multiqc_report.html",
        config['paths']['output'] + "/qc_report/qc-report.html",
        config['paths']['output'] + "/test_associations/summary-stats"
        
