# Run analysis of a small number of TRs

To validate the wdl:

```
# Main WDL for targeted TR genotyping
java -jar womtool-84.jar validate wdl/workflows/targetTR.wdl

# Helper WDL for merging mega-batches from UKB
java -jar womtool-84.jar validate wdl/workflows/merge_index_str.wdl

# WDL for individual association testing
java -jar womtool-84.jar validate wdl/tasks/associatr.wdl
```

To test locally:

```
# One batch of 3 samples
java -jar cromwell-84.jar run -i tests/test_local.json wdl/workflows/targetTR.wdl

# 3 batches of 1 sample each
java -jar cromwell-84.jar run -i tests/test_local_2.json wdl/workflows/targetTR.wdl

# Test association tests+vis
java -jar cromwell-84.jar run -i tests/test_associatr.json wdl/tasks/associatr.wdl
```