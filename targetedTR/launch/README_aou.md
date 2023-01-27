# TargetTR on All of Us workbench

All of this must be run on the workbench itself. Once in the workbench, launch a terminal session and to the following.

1. Navigate to the `launch` directory.

2. Do `git pull` to get the latest changes.

3. Copy the `manifest.csv` file with the CRAM paths to the launch directory.
```
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v6/wgs/cram/manifest.csv .
```

4. Run a test job:

```
./targetTR_launcher_aou.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 1 \
  --batch-num 2 \
  --file-list manifest.csv
```