# Extracting a basic dataset from [MethylKit](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r87) output files

In cases where reads have already been aligned using another pipeline
and DNA methylation levels have been extracted using
[MethylKit](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r87),
we provide a script
[scripts/methylkit-to-matrix.sh](scripts/methylkit-to-matrix.sh)
to merge these outputs into an analyzable dataset that
includes the following output files:

* `cell-counts.csv` Cell count estimates per sample per cell type (assumes that the DNA methylation was measured in peripheral blood)
* `coverage.csv` Number of reads containing each CpG site for each sample
* `illumina.csv` Methylation levels per sample per CpG site contained on the Illumina 450K or EPIC beadchip
* `methylation-scores.csv` A matrix of published methylation scores per sample calculated from the methylation levels (see [meffonym](https://github.com/perishky/meffonym) for details)
* `methylation.csv` A matrix of methylation levels per sample per CpG site for all sites contained in at least 10 reads for at least 50% of the samples
* `qc.html` A report providing an overview of the extracted data

The script can be executed as follows:

```
REPO_DIR=path/to/this/repository
DATA_DIR=path/to/data/directory
PANEL=path/to/panel
OUT_DIR=path/to/output/directory
ASSEMBLY=hg38

IMAGE_URL=oras://docker.io/onuroztornaci/methylhead-pipeline:meth_analysis
IMAGE=methylhead-pipeline_meth_analysis.sif
if [ ! -e $IMAGE ]; then
    apptainer pull $IMAGE_URL
fi
apptainer exec -B $REPO_DIR -B $DATA_DIR -B $OUT_DIR $IMAGE \
    $REPO_DIR/scripts/methylkit-to-matrix/methylkit-to-matrix.sh \
    $REPO_DIR \
    $DATA_DIR \
    $PANEL \
    $OUT_DIR \
    $ASSEMBLY
```

The 'panel' defines the genomic target regions formatted as a BED file
(i.e. at least tab-separated columns corresponding to chromosome, start and end).
For example, the Twist methylome panel BED files can be found here:
https://www.twistbioscience.com/products/ngs/fixed-panels/human-methylome-panel?tab=resources

The final argument for the script, the genome assembly is optional.
The default is 'hg19'.  If read alignment was to the hg38 assembly, then
'hg38' needs to be specified.

> The script 'methylkit-to-matrix.sh' can be run outside of a container. 
> It requires that R is installed with the following R packages:
> * IlluminaHumanMethylation450kanno.ilmn12.hg19
> * IlluminaHumanMethylationEPICanno.ilm10b4.hg19
> * IlluminaHumanMethylationEPICv2anno.20a1.hg38
> * data.table
> * dplyr
> * methylKit
> * [meffonym](https://github.com/perishky/meffonym)
