# cast-workflows
Workflows used by the Center for Admixture Science and Technology

To compile targetedTR for DNANexus, run from the current directory:

```
java -jar "$DX_COMPILER_JAR" compile \
        targetedTR/wdl/workflows/targetTR.wdl \
        -project project-GG25fB8Jv7B928vqK7k6vYY6 \
        -folder "$PLACE_TO_PUT_COMPILED_APP" \
        -streamFiles all \
        -archive
```
