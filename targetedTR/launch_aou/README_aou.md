# TargetTR on All of Us workbench

All of this must be run on the workbench itself. Once in the workbench,start both cromwell and cloud environment. Open notebook and add cromwell setup snippet. Run the cell and launch a terminal session and do the following.
 
1. Clone the github and navigate to the launch directory.

```
cd workspaces/impactofglobalandlocalancestryongenomewideassociationv6studies/
mkdir mgymrek-workspace
cd mgymrek-workspace/
git clone https://github.com/cast-genomics/cast-workflows/
cd cast-workflows/
git checkout nichole/cromshell
cd targetedTR/launch_aou
```

2. Set up (already set up in jupyter snippet)

3. Copy the `manifest.csv` file with the CRAM paths to the launch directory.
```
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v6/wgs/cram/manifest.csv .
```

4. Run a test job:

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
