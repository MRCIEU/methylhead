import os
import subprocess
import pandas
from snakemake.logging import logger

# ===================================================================
# File paths
# ===================================================================

# === normalize_path(path, base): return absolute path with no shortcuts 
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

# === dirname(path): return directory of file path =================
def dirname(path):
    return os.path.dirname(normalize_path(path))

# === normalize config file paths ===========================
for name,path in config['inputs'].items():
    config['inputs'][name] = normalize_path(path)

for name,path in config['paths'].items():
    config['paths'][name] = normalize_path(path)
    
# === config:inputs:fasta--path to genome assembly fasta file ==============
config['inputs']['fasta'] = f"{config['paths']['genome']}/{config['parameters']['assembly']}.fa"

# === config:paths:scripts--directory of repository scripts ===================
config['paths']['scripts'] = workflow.basedir + "/scripts"
    
# === create output directory =========================================
try:
    os.makedirs(config['paths']['output'])
except FileExistsError:
    pass
    
# =====================================================================
# Samples
# =====================================================================

# === samples: samplesheet ============ ===============================
samples = pandas.read_csv(config["inputs"]["samplesheet"])

# === all_samples: list of sample ids ===============================
all_samples = samples["sample_id"].unique().tolist()

# === normalize fastq file paths =====================================
fastq_dir = config["inputs"]["fastq"]
samples["read1"] = [normalize_path(path, fastq_dir) for path in samples["read1"]]
samples["read2"] = [normalize_path(path, fastq_dir) for path in samples["read2"]]

# =====================================================================
# Apptainer 
# =====================================================================

# === pull container images ==========================================
def pull_image(url,image_file):
    if not os.path.exists(image_file):        
        subprocess.run(["apptainer","pull",image_file,url],check=True)
    return image_file
for container in config['containers'].values():
    pull_image(container['image_url'], container['image_file'])
    
# === apptainer_exec(container,cmd): generate command to run 'cmd' in container ===== 
def apptainer_exec(container, cmd, extra_args="--writable-tmpfs"):
    bind_paths = [path for path in config['paths'].values()] \
        + [dirname(path) for path in config['inputs'].values()] \
        + [os.getcwd()]
    bind_paths = list(set(bind_paths))
    apptainer_args = (
        " ".join([ f"-B {path}:{path}" for path in bind_paths ]) +
        " --pwd " + os.getcwd())
    image = config['containers'][container]['image_file']
    return f"apptainer exec {apptainer_args} {extra_args} {image} /bin/bash -c \"{cmd}\""

# =====================================================================
# HPC resources
# =====================================================================

# === get_resources(): hpc resources required =============================
def get_resources(resource_class="light", **overrides):
    resources = config["resources"][resource_class].copy()
    resources.update(overrides)
    return resources
        
# ======================================================================
# CHECKS
# ======================================================================

# === check config['inputs']['samplesheet'] is csv with columns sample_id,read1,read2 =
# === check config['inputs']['phenotypes'] is csv with column sample_ id ==============
# === check config['inputs']['models'] is a csv file with columns name,var,model ======
# === check config['inputs']['panel'] is csv with columns chr,start,end ===============
# === check config['inputs']['region_model'] is text file with a single line defining a model with variables in pheno =====

# === check fastq files exist and are not empty =============================
def is_fastq_ok(filename, min_size=2**20):
    if not os.path.exists(filename):
        logger.error(f"File not found: {filename}")
        return False
    size = os.path.getsize(filename)
    if size < min_size:
        logger.info(f"File is too small: {filename}\n  {round(size/2**20,1)}Mb")
        return False
    return True
fastq_files = samples['read1'].to_list() + samples['read2'].to_list()
if not all([is_fastq_ok(fq) for fq in fastq_files]):
    print("There are problems with one or more fastq files, see details above.\n")
    sys.exit(1)

# ===================================================================
# Rule modifications
# ===================================================================

test_associations = 'models' in config['inputs']

identify_regions = 'region_model' in config['inputs']

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
include: "rules/methyldackel_bedgraph.smk"
include: "rules/meth_bedgraph.smk"
include: "rules/extract_methyl_entropy.smk"
include: "rules/extract_methylkit.smk"
include: "rules/methylation_matrix.smk"
include: "rules/illumina_matrix.smk"
include: "rules/methyl_entropy_matrix.smk"
include: "rules/estimate_cell_counts.smk"
include: "rules/methylation_scores.smk"
include: "rules/raw_read_counts.smk"
include: "rules/qc_report.smk"

if identify_regions:
   include: "rules/methylseqlm.smk"

if test_associations:
   include: "rules/test_associations.smk"

rule all:
    input:
        config['paths']['output'] + "/methylation/matrix.csv.gz",
        config['paths']['output'] + "/illumina/matrix.csv.gz",
        config['paths']['output'] + "/methyl_entropy/matrix.csv.gz",
        config['paths']['output'] + "/cell_counts/matrix.csv",
        config['paths']['output'] + "/methylation_scores/matrix.csv",
        config['paths']['output'] + "/multiqc/multiqc_report.html",
        config['paths']['output'] + "/qc_report/qc-report.html",
        config['paths']['output'] + "/methylseqlm/matrix.csv.gz" if identify_regions else [],
        config['paths']['output'] + "/test_associations/summary-stats" if test_associations else []
        
        
        
