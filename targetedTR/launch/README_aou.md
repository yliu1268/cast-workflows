# TargetTR on All of Us workbench

All of this must be run on the workbench itself. Once in the workbench, launch a terminal session and to the following.

1. Navigate to the `launch` directory.

2. Do `git pull` to get the latest changes.

3. Copy the `manifest.csv` file with the CRAM paths to the launch directory.

4. Run a test job:

```
./targetTR_launcher_aou.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 5 \
  --batch-num 3 \
  --file-list manifest.csv
```