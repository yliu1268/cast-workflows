# TargetTR on UKB RAP

## Setup WDL workflow

First, compile WDL workflow to run on DNA Nexus

```
# Notes
# - Checks are more rigorous than those of womtool validate
java -jar ~/Applications/dxCompiler-2.10.4.jar compile ../wdl/workflows/targetTR.wdl \
	-project project-GG25fB8Jv7B928vqK7k6vYY6 -folder /TargetedSTR/ -streamFiles all -archive 
```

## Compile list of files to process

List all UKB cram/crai files
```
for i in $(seq 10 60)
do
  cmd="dx ls -l 'Bulk/Whole genome sequences/Whole genome CRAM files/${i}/'"
  sh -c ${cmd}
done > ukb_cram_files_long.txt

./process_cram_list.py ukb_cram_files_long.txt > ukb_cram_and_index_files.txt
```

## Call python launcher 

The python launcher takes care of batching, which was difficult to do cleanly in WDL 1.0.

The code below shows a small example (3 batches of 25 samples each).
To run all files, increase `--batch-size` to a larger number (1000?) and remove the `--batch-num` argument.
```
./targetTR_launcher_ukb.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 25 \
  --batch-num 3 \
  --workflow-id workflow-GPYY0p0Jv7B730qb33V6jvG0 \
  --file-list ukb_cram_and_index_files.txt
```

See `launch_scripts/` for full launches.