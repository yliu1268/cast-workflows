# Debug mergeSTR step

WORK IN PROGRESS - instructions to run only mergeSTR step from a prevous targetTR attempt

```
name="" # set to job name from targetTR --name option

# Set up input json
# needs merge_hipstr.vcfs, merge_hipstr.vcf_indices, out_prefix


# Run mergeSTR step
cromshell submit ../wdl/merge_hipstr.wdl \
  debug-${name}.json \
  -d "${name}-wdl.zip" 
```