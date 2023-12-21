# Run analysis of a small number of TRs

To validate the wdl:

```
# Main WDL for targeted TR genotyping
java -jar womtool-84.jar validate wdl/targetTR.wdl

# Helper WDL for merging mega-batches from UKB
java -jar womtool-84.jar validate wdl/merge_index_str.wdl
```

To test locally:

```
# One batch of 3 samples
java -jar cromwell-84.jar run -i tests/test_local.json wdl/targetTR.wdl

# 3 batches of 1 sample each
java -jar cromwell-84.jar run -i tests/test_local_2.json wdl/targetTR.wdl

# 3 batches of 1 sample each - string input as in AOU analyses
# This isn't working... file paths change in the dockers
java -jar cromwell-84.jar run -i tests/test_local_strings.json wdl/targetTR.wdl
```