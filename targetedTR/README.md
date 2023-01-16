# Run analysis of a single TR

To validate the wdl:

```
java -jar womtool-84.jar validate wdl/workflows/targetTR.wdl
```

To test locally:

```
java -jar cromwell-84.jar run -i tests/test_local.json wdl/workflows/targetTR.wdl
```