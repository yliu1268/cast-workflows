**UKB launcher under construction**

# TODO

* where are the reads in the new release??
* check final batch error
* Make take a bed file rather than region, similar to AoU  launcher
* use new wdl arguments

# TargetTR on UKB RAP

The steps below can be run from the command line on your local computer (or anywhere else you prefer to run and test code) to launch jobs on DNA Nexus. You can use the DNA Nexus web console to track job status and browse files.

(Note, this is in contrast to All of Us where all steps must be run on the AoU workbench).

## Prerequisites

* Make sure the [DNA Nexus toolkit](https://documentation.dnanexus.com/downloads) is installed. You can install it with `pip3 install dxpy`
* Run `dx login` to log in with your DNA Nexus credentials.
* Download the jar file for [dxCompiler](https://github.com/dnanexus/dxCompiler/releases). This is a tool that compiles WDLs to convert them into DNA Nexus-compatible workflows.

## Setup WDL workflow

Running targetTR on DNANexus uses two workflows:

* The main targetTR workflow defined in `../wdl/targetTR.wdl`
* An auxiliary workflow defined in `../wdl/merge_index_str.wdl`

This is because launching a single targetTR run on all files at once stalls on the command line, probably because the JSON to describe the job is too big. To get around this, the launcher script runs multiple targetTR jobs then does a final step to merge the results.

First, compile WDL workflows to run on DNA Nexus

```
DX_COMPILER_JAR=dxCompiler-2.10.4.jar # change to path to your jar file
PROJID=project-GG25fB8Jv7B928vqK7k6vYY6 # change to appropriate project ID if needed
# Note: Checks are more rigorous than those of womtool validate
java -jar "$DX_COMPILER_JAR" compile ../wdl/targetTR.wdl \
	-project $PROJID -folder /TargetTR/ -streamFiles all -archive -separateOutputs

java -jar "$DX_COMPILER_JAR" compile ../wdl/merge_index_str.wdl \
  -project $PROJID -folder /TargetTR/ -streamFiles all -archive 
```

## Compile list of files to process

The steps below need to run once to get a list of samples to process. These steps do the following:

* List all UKB cram/crai files
* Remove problematic samples (**Note! this is important. some samples seem to have correct CRAM files. If included they will break the whole job.** These are removed by `process_cram_list.py`)
* Process that into a file that contains the IDs of the CRAM and index files for each sample to be processed.

```
for i in $(seq 10 60)
do
  cmd="dx ls -l 'Bulk/Whole genome sequences/Whole genome CRAM files/${i}/'"
  sh -c "${cmd}"
done > ukb_cram_files_long.txt

./process_cram_list.py ukb_cram_files_long.txt > ukb_cram_and_index_files.txt
```

## Run a small test job

The code below shows a small example (3 batches of 25 samples each).
```
./targetTR_launcher_ukb.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 25 \
  --batch-num 3 \
  --workflow-id workflow-GPfbXV8Jv7B27kpf6Y50QyQ9 \
  --file-list ukb_cram_and_index_files.txt
```

## Run a full job on all samples
TODO

See `launch_scripts/` for full launches.

## Detailed usage for targetTR_launcher_ukb.py

Required options:
* `--name <STRING>`: Name of the run. Used as the prefix to output files.
* `--tr-bed <PATH>`: Path (local file or DNA Nexus file ID) to BED file of TRs to run (HipSTR reference format). 
* `--file-list <PATH>`: Path (local file) with list of samples to be processed. Format of each line: cram-file-id cram-index-id. This is created above by `process_cram_list.py`.

Additional input files:
* `--genome-id <PATH>`: DNANexus file ID of reference genome. Defaults to a path to a version of `Homo_sapiens_assembly38.fasta`.
* `--genome-idx-id <PATH>`: DNANexus file ID of the index (`.fai`) file of the reference genome. Defaults to a version of `Homo_sapiens_assembly38.fasta.fai`

Additional run options:
* `--workflow-id <STR>`: ID of the DNANexus workflow for targetTR. This is set to a default value which was the working version of the workflow at the time of writing. **If you modify the workflow, be sure to change the workflow ID to make sure the new one is being used!**.
* `--batch-size <INT>`: Number of samples to process in each batch. Default: 300.
* `--batch-num <INT>`: Number of batches to process. Default: -1 (process all batches). This is helpful to set during debugging to consider a small number of batches.
* `--merge-workflow-id <STR>`: DNANexus workflow ID of the merge index workflow. Defaults to a recent working version.
* `--max-batches-per-workflow <INT>`: Maximum number of batches to launch at once. Default: -1 (all)
* `--concurrent`: Launch all batches at once. (otherwise, launch one at a time)

