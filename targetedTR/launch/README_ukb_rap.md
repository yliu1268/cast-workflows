Compile to run on DNA Nexus

```
# Notes
# - Checks are more rigorous than those of womtool validate
java -jar ~/Applications/dxCompiler-2.10.4.jar compile ../wdl/workflows/targetTR.wdl \
	-project project-GG25fB8Jv7B928vqK7k6vYY6 -folder /TargetedSTR/ -streamFiles all -archive 
```

Run example

```
# Test one batch of 3
dx run workflow-GP8jzvQJv7B3K97ypgvBBqxq -y -f uk_rap_example_input.json --destination "TargetedSTR/results/CSTB"

# Test three batches of 1
dx run workflow-GP8jzvQJv7B3K97ypgvBBqxq -y -f uk_rap_example_input-2.json --destination "TargetedSTR/results/CSTB-2"

```

List all UKB cram files
```
for i in $(seq 10 60)
do
	cmd="dx ls 'Bulk/Whole genome sequences/Whole genome CRAM files/${i}/'"
	sh -c ${cmd} | awk -v "i=$i" '{print "Bulk/Whole genome sequences/Whole genome CRAM files/" i "/" $0}'
done | grep -v "crai$" > ukb_cram_files.txt
```

Get file IDs of all the crams/indices
TODO need to rerun this takes a long time
```
while IFS="" read -r p || [ -n "$p" ]
do
	cramid=$(dx describe "$p" | awk '($1=="ID") {print $2}')
	idxid=$(dx describe "$p.crai" | awk '($1=="ID") {print $2}')
  	echo $cramid $idxid
done < ukb_cram_files.txt > ukb_cram_and_index_files.txt
```

Run bigger example (3 batches of 100 each)
```
./targetTR_launcher_ukb.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 5 \
  --batch-num 3 \
  --workflow-id workflow-GP8jzvQJv7B3K97ypgvBBqxq \
  --file-list ukb_cram_and_index_files.txt
```
