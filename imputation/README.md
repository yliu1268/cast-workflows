# Run imputation on AoU dataset

The goal of imputation is to perform phasing and imputation on a cohort of samples. The workflow performs the following steps:

1. **Genotype one or more STRs using HipSTR.** TargetTR performs multi-sample calling in batches. This is because using HipSTR to jointly genotype a huge number of samples (e.g. 200,000) is too memory intensive. After genotyping each batch, it uses mergeSTR to combine the output VCF files.
2. **Perform basic quality filters with dumpSTR.** DumpSTR is invoked with the following best practice filters for HipSTR genotyping: `--hipstr-min-call-Q 0.9 --hipstr-min-call-DP 10 --hipstr-max-call-DP 10000 --hipstr-min-supp-reads 2 --hipstr-max-call-stutter 0.15 --hipstr-max-call-flank-indel 0.15`
3. **Sort, bgzip and tabix index the output VCF file.**


**For most of our use cases, users do not interact with the WDL described here directly, but rather call launcher scripts ** which handle setting up inputs and calling the WDL on particular cloud systems. Sections below give additional WDL details which can be helpful for development/testing or debugging.

