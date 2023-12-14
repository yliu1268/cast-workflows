# TargetTR on All of Us workbench

All of this must be run on the workbench itself. Once in the workbench, before running the commands below:

* Start the cromwell environment (pink circle icon on the right).
* Start a cloud environment (Jupyter icon on the right). The instructions below were tested tih the "Recommended environment" with 4 CPUs/26GB RAM.
* Navigate to the "Analysis" tab in the workspace and open the "cromshell.ipynb" notebook. Run the first cell to configure cromshell. 
* Open a terminal (terminal icon on the right) and run the commands below.
 
1. Clone the github and navigate to the launch directory.

```
cd workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies
mkdir mgymrek-workspace
cd mgymrek-workspace/
git clone https://github.com/cast-genomics/cast-workflows/
cd cast-workflows/
git checkout mgymrek/cromshell3
cd targetedTR/launch_aou
```

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

Next attempts:
* print out chip/hardware info being used to see if it correlates with segfaults?
* catch 139 (segfault) gracefully and still create good (empty) vcf+index

Notes on running full:
*  seems to get a loy of "RetryableFailure" with segfault that then succeed on later tries. only happens when we run a huge amount at once. wonder if we can decrease batch size, or CPU?, or number of jobs run simultaneously? if that would help. maybe we overwhelm the system when we submit >800 jobs in parallel?
* what is memoryMin and cpuMin in WDL runtime?

Weird stuff for cromshell logs:
* "No logs with status ['ALL'] found for workflow, try adding the argument '-s ALL' to list logs with any status"
* cromshell list status change after I run cromshell status on a particulra job id
* cromshell logs -s ALL seem to be missing at least one line each time in output

After this works on AoU:
* rename aou/ukb options to be more general
* test new aou workflow on ukb
* document all wdl options
* merge to master branch of cast-workflows
* github actions for dockers? also move dockers up one directory level so they can be shared by different pipelines?
