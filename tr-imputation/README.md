# Running TR imputation on the All of Us workbench

The goal of imputation is to perform phasing and imputation of TRs on a cohort of samples. The workflow performs the following steps:

1. **TR imputation with Beagle.** We perform phasing with Beagle in batches. This is because using running Beagle on a large number of samples (e.g. 200,000) is too memory intensive.
2. **Merge imputed TR and SNPs seperately with bcftools.** 
3. **Sort, bgzip and tabix index the output VCF file.**
4. **Run annotaTR to add dosages and convert to PGEN format.**

For most of our use cases, users do not interact with the WDL described here directly, but rather call launcher scripts which handle setting up inputs and calling the WDL using cromshell. Sections below give additional WDL details which can be helpful for development/testing or debugging.

## Setup
In all cases you'll have to run the following steps in the AoU workbench before starting a workflow:

1. Start the cromwell environment (pink circle icon on the right).
2. Start a cloud environment (Jupyter icon) with:
    * "General Analysis" environment
    * 2 CPU/7.5GB RAM
    * Compute type: Standard VM
    * Automatically pause after: 30 minutes
    * Storage disk options: standard disk, 120GB
3. Open a terminal (terminal icon) and run the commands below:

```
git clone https://github.com/cast-genomics/cast-workflows/
cd cast-workflows/tr-imputation
../utils/configure-cromshell.py
```

## Run a small test job

It is recommended to first run a small test on a couple samples to make sure everything is set up correctly. e.g.:

```
chrom=21
./batch_imputation_aou.py \
--name chr${chrom}_test \
--batch-num 2 \
--vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
--chrom ${chrom}

```

This will print out a friendly turtle with the job ID if successful. Use the following to check the status of your job. It will take around 10-20 minutes to run. If successful the status will eventually change to "Succeeded".

```
cromshell status $JOBID
```

If you see "Failed", you can look at the logs to see what happened:

```
cromshell logs -s ALL $JOBID
```

You can check the output:
```
cromshell list-outputs $JOBID
```

## Run a full job on all samples

```
chrom=21
./batch_imputation_aou.py \
--name chr${chrom} \
--vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
--chrom ${chrom}
```

To extract the Beagle files and upload separately to `${WORKSPACE_BUCKET}/beagle_hg19/chr${chrom}`

```
cromshell -mc list-outputs -j -d $jobid | python -c "import json, sys; data=json.load(sys.stdin); [sys.stdout.write(item['outvcf']+'\n'+item['outvcf_index']+'\n') for item in data['batch_imputation.beagle']]" | xargs -n1 -I% -P1 sh -c "gsutil mv % ${WORKSPACE_BUCKET}/beagle_hg38/chr${chrom}"
```

## To check on status on your run
```
cromshell logs -s ALL -des $JOBID (use -des for descriptive log info)
cromshell logs -s Failed -des $JOBID
```
