Compile to run on DNA Nexus

```
# Notes
# - Checks are more rigorous than those of womtool validate
java -jar ~/Applications/dxCompiler-2.10.4.jar compile ../wdl/workflows/targetTR.wdl \
	-project project-GG25fB8Jv7B928vqK7k6vYY6 -folder /TargetedSTR/ -streamFiles all -archive
```

Run example

```
dx run workflow-GP3JQGQJv7B2q9PG5fk6Kqg8 -y -f uk_rap_example_input.json --destination "TargetedSTR/results/CSTB"
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
```
dx describe file-G2x1p68JkF61f5xgGk9XyffQ | awk '($1=="ID") {print $2}'
```