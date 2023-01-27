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

List all UKB cram files
```
for i in $(seq 10 60)
do
	cmd="dx ls 'Bulk/Whole genome sequences/Whole genome CRAM files/${i}/'"
	sh -c ${cmd} | awk -v "i=$i" '{print "Bulk/Whole genome sequences/Whole genome CRAM files/" i "/" $0}'
done | grep -v "crai$" > ukb_cram_files.txt
```

Get file IDs of all the crams/indices
(this takes a long time. is there a faster way?)
```
while IFS="" read -r p || [ -n "$p" ]
do
	cramid=$(dx describe "$p" | awk '($1=="ID") {print $2}')
	idxid=$(dx describe "$p.crai" | awk '($1=="ID") {print $2}')
  	echo $cramid $idxid
done < ukb_cram_files.txt > ukb_cram_and_index_files.txt
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
  --workflow-id workflow-GP9f1pjJv7B4G6fk11V1y5QF \
  --file-list ukb_cram_and_index_files.txt
```
