# Nextflow for rMATS

* [rmatsTurbo.nf](rmatsTurbo.nf) is a Nextflow ([https://www.nextflow.io](https://www.nextflow.io)) implementation of the rMATS workflow with BAM files as input


## Run with nextflow

> It requires `Java` environment.

### Installation of Nextflow
* [https://www.nextflow.io/docs/latest/install.html](https://www.nextflow.io/docs/latest/install.html)

###  Run the workflow

Modify the `nextflow.config` file to set the parameters for the workflow. See the example in the `example.config` file.

```bash
nextflow run rmatsTurbo.nf --config nextflow.config -resume
```

### Run on the cloud
Nextflow was designed to run on the cloud. You can use the `nextflow.config` file to set the parameters for the cloud execution. More details can be found in the Nextflow [documentation](https://www.nextflow.io/docs/latest/aws.html).


