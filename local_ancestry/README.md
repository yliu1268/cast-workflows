# Setup

This must be run on the All of Us workbench. You must have launched both an analysis and cromshell environment. This also requires having run the `subset_vcf` workflow to generate batch VCFs at `${WORKSPACE_BUCKET}/acaf_batches/chr${chrom}`.

```
git clone https://github.com/cast-genomics/cast-workflows
cd cast-workflows/
git checkout gnomix
cd local_ancestry
../utils/configure-cromshell.py
```

# Test

```
chrom=11; ./gnomix_launcher.py \
  --name test-chr${chrom} \
  --vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
  --chrom ${chrom} \
  --max-batches 2 
```

# Full chromosome

```
chrom=11; ./gnomix_launcher.py \
  --name gnomix-chr${chrom} \
  --vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
  --chrom ${chrom}
```

# Upload outputs to our workspace bucket

```
# Get msp/fb files and copy to ${WORKSPACE_BUCKET}/gnomix/outputs/
cromshell list-outputs $jobid

# Also keep track of useful phased Beagle files!
cromshell -t 2000 -mc list-outputs -j -d $jobid > out.json 
cat out.json | python -c "import json, sys; data=json.load(sys.stdin); [sys.stdout.write(item['run_gnomix.beagle'][0]['outvcf']+'\n'+item['run_gnomix.beagle'][0]['outvcf_index']+'\n') for item in data['local_ancestry.run_gnomix']]" | xargs -n1 -I% -P1 sh -c "gsutil mv % ${WORKSPACE_BUCKET}/beagle_hg19/chr${chrom}/"
```

Old:

```
# Note: now metadata is printing the colors even when using -mc?
# use list-outputs instesad see above
chrom=chr11
cromshell -mc -t 2000 metadata $jobid > gnomix_${chrom}_metadata.json
python -c "import json; data = json.load(open('gnomix_chr${chrom}_metadata.json', 'r')); [print(data['calls']['local_ancestry.run_gnomix'][i]['subWorkflowMetadata']['calls']['run_gnomix.beagle'][0]['outputs']['outvcf']) for i in range(len(data['calls']['local_ancestry.run_gnomix']))]; [print(data['calls']['local_ancestry.run_gnomix'][i]['subWorkflowMetadata']['calls']['run_gnomix.beagle'][0]['outputs']['outvcf_index']) for i in range(len(data['calls']['local_ancestry.run_gnomix']))];" | xargs -n 1 -P 1 -I% sh -c "echo gsutil mv % ${WORKSPACE_BUCKET}/beagle_hg19/${chrom}"
```