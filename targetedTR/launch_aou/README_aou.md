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
git checkout nichole/cromshell3
cd targetedTR/launch_aou
```

2. Copy the `manifest.csv` file with the CRAM paths to the launch directory.
```
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v6/wgs/cram/manifest.csv .
```

3. Run a test job:

```
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
export GCS_REQUESTER_PAYS_PROJECT=${GOOGLE_PROJECT};

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

To check the status of your job, you can run:
```
cromshell status $JOBID
```  
       
To check medata and log, you can run:
```
cromshell -t 20 metadata $JOBID
cromshell logs -s ALL $JOBID
```
List all output files produced by a workflow, you can run:
```
cromshell list-outputs $JOBID
```
# TDL


* remove raw urls to wdl dependencies
* rename aou/ukb options to be more general
* test new aou workflow on ukb
* document all wdl options
* merge to master branch of cast-workflows
