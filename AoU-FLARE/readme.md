# Local Ancestry Calling with FLARE (chr1-22, chrX)

**Last updated**: Feb 17, 2025  
**FLARE version**: 0.5.2  

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Workflow Description](#workflow-description)
4. [Input & Output Descriptions](#input--output-descriptions)
5. [Recommended AoU Environment Settings](#recommended-aou-environment-settings)
6. [Example Usage for chr2](#example-usage-for-chr2)
7. [Example Results](#example-results)
8. [Running for chr1-22, chrX](#running-for-chr1-22-chrx)
9. [References & Citation](#references--citation)

---

## Overview

This workflow uses [**FLARE**](https://github.com/browning-lab/flare) to perform **local ancestry** calls on **phased VCF** files. It is configured for the [All of Us (AoU) Researcher Workbench](https://www.researchallofus.org/for-researchers/workbench/), where large-scale WGS data can be processed.

We provide a **WDL-based** pipeline for running FLARE efficiently. However, for large chromosomes, we have also tested an alternative **simple shell script approach**, which works well:

- **`run-flare-chr2-full-code.sh`**  
  [Sample Run] Run FLARE on Chr2 with a simple script.
- **`run-FLARE.sh`**  
  [Full Code] Run FLARE on a specific chromosome with a user-specified chunk size.  
  This script can be run directly in the **AoU cloud environment**, and with interactive prompts, users only need to **input the chromosome number** and **chunk size** to execute FLARE.

The main steps are:

1. **Split** a chromosome-wide, **BEAGLE-phased VCF** into **sample chunks** (using `bcftools`).
2. **Run FLARE** on each chunk in parallel.
3. **Merge** the local ancestry VCF files and global ancestry text files across all chunks into two final files per chromosome.

By default, the WDL **scatters** over sample chunks for a single chromosome. You can run the workflow on **chr1–22** plus **chrX** sequentially or modify the WDL to handle multiple chromosomes in one submission.

---

## Requirements

1. **Phased VCF** for each chromosome, indexed (`.vcf.gz` + `.tbi`):
   - Must be **phased** with **no missing** genotypes.  
   - Tools like [BEAGLE](https://faculty.washington.edu/browning/beagle/) can generate these.
2. **Reference panel** for local ancestry, e.g. 1000 Genomes:
   - A reference VCF file (`ref=...`).  
   - A panel file (`ref-panel=...`) mapping reference-sample IDs to population labels.
3. **Genetic map** in PLINK format (e.g., `plink.chr2.GRCh38.map`) with cM units.
4. **Java** >= 1.8 and **FLARE** jar (≥ 0.5.2).
5. **bcftools** for subsetting and merging.
6. AoU environment with **enough memory** and **CPU** to handle chunk splitting, running FLARE, and merging results (see [Recommended AoU Environment Settings](#recommended-aou-environment-settings)).
7. A **WDL** execution engine, e.g. [Cromwell](https://github.com/broadinstitute/cromwell) or [cromshell](https://github.com/broadinstitute/cromshell).

---

## Workflow Description

A typical WDL file (e.g., `local_ancestry.wdl`) defines three main tasks:

1. **ChunkVCF**  
   - Splits the input VCF by sample subsets (defined by `chunk_size`).
   - For each subset, writes `chunk_#.vcf.gz` and an index.
2. **RunFlare**  
   - Runs `java -jar flare.jar` for each chunked VCF:
     - `ref`, `ref-panel`, `gt`, `map` are required parameters.
     - `probs=true/false` to include posterior probabilities.
     - `em=true/false` to enable or disable iterative parameter estimation.
   - Produces `.anc.vcf.gz` (local ancestry) and `.global.anc.gz` (global ancestry).
3. **IndexMergeOutputs**  
   - Indexes each `.anc.vcf.gz` if needed.
   - Merges them with `bcftools merge` → produces a single `merged.anc.vcf.gz`.
   - Concatenates all `.global.anc.gz` files via `zcat` → produces one `merged.global.anc.gz`.

Final results are:

- **`<prefix>.merged.anc.vcf.gz`**  
- **`<prefix>.merged.global.anc.gz`**

---

## Recommended AoU Environment Settings

- **Chunk Splitting** & **Merging** with `bcftools`:
  - Suggest **32 CPUs** and **120 GB RAM**.
- **Running FLARE**:
  - Suggest **96 CPUs** and **624 GB RAM**.
  - FLARE can be memory-intensive, especially when generating posterior probabilities (`probs=true`).

---

## Example Usage for chr2

### 1. `local_ancestry.json`

```json
{
  "LocalAncestryWorkflow.chunk_size": 12000,
  "LocalAncestryWorkflow.chromosome": "2"
}
```

### 2. Example Run Command

```bash
cromshell submit local_ancestry.wdl local_ancestry.json
```

Alternatively, for a direct AoU environment run, use:

```bash
bash run-FLARE.sh
```

---

## Example Results

```
#CHROM   POS     ID      REF  ALT  FORMAT                 Sample1                                     Sample2
19       10000   .       A    C    GT:AN1:AN2:ANP1:ANP2   0|0:2:1:0.39,0,0.61,0,0:0.19,0.46,0.35,0,0    0|0:3:0:0.38,0,0.01,0.62,0:0.99,0,0,0,0
```

---

## References & Citation

- **FLARE** developed by the Browning lab:  
  [https://github.com/browning-lab/flare](https://github.com/browning-lab/flare)
