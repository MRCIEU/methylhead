**Picard Usage:**

```
nextflow Picard_pipeline.nf \
  --reads "[fastq path]/*_R{1,2}*.fastq.gz" \
  --intervals covered_targets_Twist_Methylome_hg19_annotated_collapsed_final
  --genome_folder  [BWA genome index path]
 preview -with-dag flowchart.html \
  -resume \
  -with-timeline time_line.html \
  -with-report report.html

```
preview -with-dag flowchart
 
 The pipeline will be represented as a direct acyclic graph (DAG)

- with-timeline
 
 Using to enable the creation of the timeline report.

- with-report
 
 It creates an HTML execution report: a single document about resources usage (which includes many useful metrics about a workflow execution).
