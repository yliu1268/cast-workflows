# Run analysis of a single TR

To validate the wdl:

```
java -jar womtool-84.jar validate wdl/workflows/targetTR.wdl
java -jar womtool-84.jar validate wdl/workflows/merge_index_str.wdl

```

To test locally:

```
# One batch of 3 samples
java -jar cromwell-84.jar run -i tests/test_local.json wdl/workflows/targetTR.wdl

# 3 batches of 1 sample each
java -jar cromwell-84.jar run -i tests/test_local_2.json wdl/workflows/targetTR.wdl

```