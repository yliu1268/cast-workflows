# Targeted genotyping of a small number of STRs in large cohorts

The goal of TargetTR is to perform tareted genotyping of one or more TRs from next-generation sequencing data for a cohort of samples. The workflow performs the following steps:

1. **Genotype one or more STRs using HipSTR.** TargetTR performs multi-sample calling in batches. This is because using HipSTR to jointly genotype a huge number of samples (e.g. 200,000) is too memory intensive. After genotyping each batch, it uses mergeSTR to combine the output VCF files.
2. **Perform basic quality filters with dumpSTR.** DumpSTR is invoked with the following best practice filters for HipSTR genotyping: `--hipstr-min-call-Q 0.9 --hipstr-min-call-DP 10 --hipstr-max-call-DP 10000 --hipstr-min-supp-reads 2 --hipstr-max-call-stutter 0.15 --hipstr-max-call-flank-indel 0.15`
3. **Sort, bgzip and tabix index the output VCF file.**

**For most of our use cases, users do not interact with the WDL described here directly, but rather call launcher scripts (see the [launch_ukb](launch_ukb/README.md) and [launch_aou](laungh_aou/README.md) instructions)** which handle setting up inputs and calling the WDL on particular cloud systems. Sections below give additional WDL details which can be helpful for development/testing or debugging.

## WDL Inputs

WDL inputs are specified in a JSON file. Below is an example targetTR input json. More details on options are described below:

```
{
  "targetTR.cram_file_batches": [["tests/HG00138.cstb.cram","tests/NA12878.cstb.cram","tests/HG00190.cstb.cram"]],
  "targetTR.cram_index_batches": [["tests/HG00138.cstb.cram.crai","tests/NA12878.cstb.cram.crai","tests/HG00190.cstb.cram.crai"]],
  "targetTR.genome": "tests/chr21.fa",
  "targetTR.genome_index": "tests/chr21.fa.fai",
  "targetTR.outprefix": "CSTB",
  "targetTR.tr_bed": "tests/CSTB-mini.bed"
}
```

### Specifying input sequencing files:

Users must specify input sequencing files (BAM or CRAM) to be processed in one of two ways:
* `Array[Array[File]] cram_file_batches` and `Array[Array[File]] cram_index_batches` are lists of lists of CRAM file paths and CRAM file index paths. Each internal list contains files to be processed in a single batch. e.g. to process 6 CRAM files in 2 batches of 3, you would specify:

```
cram_file_batches = [["cram1.cram", "cram2.cram", "cram3.cram"], \
	["cram4.cram", "cram5.cram", "cram6.cram"]]
cram_index_batches = [["cram1.crai", "cram2.crai", "cram3.crai"], \
	["cram4.crai", "cram5.crai", "cram6.crai"]]
```

  Note when using this set of options, by default files are copied over to the instance where the WDL is running to be processed. If using DNA Nexus, files can be streamed instead of copying. It is not recommended to use these options when streaming is not available, since it will have large and expensive space requirements.

* `Array[File] cram_file_batches_str`: provides a list of text files, where each text file contains a list of paths (e.g. google cloud paths) to CRAM or BAM files to be processed (one file per line). Notes:
  * If using this option, files are not copied over to the instance. Instead, the strings of the filenames are passed to HipSTR genotyping. TargetTR invokes a version of HipSTR that can read directly from cloud paths.

  * This option assumes index files are available at the expected cloud paths (e.g. if `gs://mysample.cram` is one of the filenames, it will expect that `gs://mysample.crai` is available).

  * if reading from non-public Google cloud buckets, you must also set the `String GOOGLE_PROJECT` and `String GCS_OAUTH_TOKEN` variables.

### Additional required inputs:

* `String outprefix`: Output files all start with this string

* `File genome`: path to reference fasta file

* `File genome_index`: path to fasta file index

* `File tr_bed`: path to bed file used for HipSTR genotyping

### Additional optional inputs:

* `Boolean infer_samps_from_file`: If set to `true`, use the BAM/CRAM filenames to infer sample names rather than relying on the read group tag. Set to true for UKB.
* `Float sleep_const`: To avoid launching too many jobs at once, HipSTR jobs can sleep for a bit before running. The number of seconds to sleep for each batch is set to `sleep_constant*batch_num`. If `sleep_const` is 0, jobs will not sleep.

## WDL Outputs 

TargetTR outputs:
* `$outprefix.filtered.vcf.gz`: Sorted and bgzipped VCF file with TR genotypes
* `$outprefix.filtered.vcf.gz.tbi` VCF index

## Testing the WDL

Note the validation steps below assume you have the files `womtool-84.jar` and `cromwell-84.jar` in the current directory.

To validate the wdl:

```
# Main WDL for targeted TR genotyping
java -jar womtool-84.jar validate wdl/targetTR.wdl
```

To test locally (note the `test/chr21.fa` file is not in the github because it is too large):

```
# One batch of 3 samples
java -jar cromwell-84.jar run -i tests/test_local.json wdl/targetTR.wdl

# 3 batches of 1 sample each
java -jar cromwell-84.jar run -i tests/test_local_2.json wdl/targetTR.wdl

# Test optional infer_samps_from_file and sleep_constant
java -jar cromwell-84.jar run -i tests/test_local_3.json wdl/targetTR.wdl

# String input as in AOU analyses
java -jar cromwell-84.jar run -i tests/test_local_strings.json wdl/targetTR.wdl
```