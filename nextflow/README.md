# Nextflow for rMATS

* [rmatsTurbo.nf](rmatsTurbo.nf) is a Nextflow ([https://www.nextflow.io](https://www.nextflow.io)) implementation of the rMATS workflow with BAM files as input

## Run

### Locally

* Install: [https://www.nextflow.io/docs/latest/install.html](https://www.nextflow.io/docs/latest/install.html)
* Edit [nextflow.config](nextflow.config)
* Run:
```bash
nextflow run rmatsTurbo.nf --config nextflow.config -resume
```

### On the cloud

Nextflow can be configured to run on the cloud. More details can be found in the Nextflow [documentation](https://www.nextflow.io/docs/latest/aws.html). Another option is to use [https://seqera.io/platform/](https://seqera.io/platform/). With a compute environment set up in Seqera, rMATS can be added as a pipeline using the URL (https://github.com/Xinglab/rmats-turbo). When launching the pipeline, most parameters are available to edit based on [../nextflow_schema.json](../nextflow_schema.json). The `bam_g1` and `bam_g2` file array parameters can be added after the other parameters:
* Click "Launch settings"
* Append to "Pipeline parameters" just before the final `}`
```
, "bam_g1": ["/path/to/file_1.bam", "/path/to/file_2.bam"], "bam_g2": ["/path/to/file_3.bam", "/path/to/file_4.bam"]
```
