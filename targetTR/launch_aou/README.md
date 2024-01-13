# Running TargetTR on the All of Us workbench

## Setup

In all cases you'll have to run the following steps in the AoU workbench before starting a workflow:

1. Start the cromwell environment (pink circle icon on the right).
2. Start a cloud environment (Jupyter icon) with:
  * "General Analysis" environment
  * 2 CPU/7.5GB RAM
  * Compute type: Standard VM
  * Automatically pause after: 30 minutes
  * Storage disk options: standard disk, 120GB
3. Open a terminal (terminal icon) and run the commands below.
4. Clone the github, navigate to the launch directory, and set up cromshell:

```
git clone https://github.com/cast-genomics/cast-workflows/
cd cast-workflows/targetTR/launch_aou
./configure-cromshell.py
```

## Run a small test job

It is recommended to first do a small test on a couple samples to make sure everything is set up correctly. e.g.:

```
./targetTR_launcher_aou.py \
  --tr-bed test.bed \
  --name mytest \
  --batch-size 2 \
  --batch-num 2 
```

This will print out a friendly turtle with the job ID if successful. Use the following to check the status of your job. It will take around 10 minutes to run. If successful the status will eventually change to "Succeeded".

```
cromshell status $JOBID
```

## Run a full job on all samples

TODO - relies on precomputed batches of size 300
to create new batch size see below

## Full usage

TODO

## Modifying the batch size

TODO

## code TDL
- default manifest file so no need to copy
- take in bed file rather than TR details
- default batch size


2. Copy the `manifest.csv` file with the CRAM paths to the launch directory.
```
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/cram/manifest.csv .
```

3. Run a test job:

```
./targetTR_launcher_aou.py \
  --region chr11:119206290-119206323 \
  --period 3 \
  --refcopies 11.0 \
  --name CBL-mini \
  --batch-size 2 \
  --batch-num 2 \
  --file-list manifest.csv 
```
(AoU use batch-size 300)

Here is a command for a full run
```
# First run the below with --action create-batches
# Then subsequent runs can skip that step with
# --action run-batches as long as 
# --batch-size and --name are the same
# TODO - change to require only batch-size is the same
# so we can reuse batches across jobs
./targetTR_launcher_aou.py \
  --region chr11:119206290-119206323 \
  --period 3 \
  --refcopies 11.0 \
  --name CBL-mini \
  --batch-size 300 \
  --file-list manifest.csv \
  --action run-batches
```

To check the status of your job, you can run:
```
cromshell status $JOBID
```  
       
To check medata and log, you can run:
```
cromshell -t 20 metadata $JOBID (increase -t timeout for large batches)
cromshell slim-metadata $JOBID
cromshell logs -s ALL $JOBID

```
List all output files produced by a workflow, you can run:
```
cromshell list-outputs $JOBID
```

To get the summarized status of all jobs in the workflow:
```
cromshell count $JOBID
```
# TDL

Weird stuff for cromshell logs:
* "No logs with status ['ALL'] found for workflow, try adding the argument '-s ALL' to list logs with any status"
* cromshell list status change after I run cromshell status on a particulra job id
* cromshell logs -s ALL seem to be missing at least one line each time in output

After this works on AoU:
* test new aou workflow on ukb
