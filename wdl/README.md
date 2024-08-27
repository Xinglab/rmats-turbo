# Workflow Description Language (WDL) for rMATS

* [rmatsTurbo.wdl](rmatsTurbo.wdl) is a WDL ([https://openwdl.org](https://openwdl.org)) implementation of the rMATS workflow with BAM files as input

## Run on AnVIL

* [https://anvil.terra.bio](https://anvil.terra.bio)
* Create a new workflow using the AnVIL website and paste or upload [rmatsTurbo.wdl](rmatsTurbo.wdl)

## Run with cromwell

* [https://github.com/broadinstitute/cromwell](https://github.com/broadinstitute/cromwell)
* Create an input config file
  + `java -jar "${WOMTOOL_JAR}" inputs --optional-inputs true rmatsTurbo.wdl > rmatsTurbo_inputs.json`
  + Edit `rmatsTurbo_inputs.json`
* Validate
  + `java -jar "${WOMTOOL_JAR}" validate --inputs rmatsTurbo_inputs.json rmatsTurbo.wdl`
* Run
  + `java -jar "${CROMWELL_JAR}" run --inputs rmatsTurbo_inputs.json rmatsTurbo.wdl`
