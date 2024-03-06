import os
import pandas as pd
import subprocess
import sys

#default params
num_batches = 5
batch_file_name = "wgs_ratios_batch"

if len(sys.argv) == 2:
    num_batches = sys.argv[1]
if len(sys.argv) == 3:
    batch_file_name = sys.argv[2]
#TODO add argument for input location

bucket = os.getenv('WORKSPACE_BUCKET')

dfs = []
for i in range(num_batches):
    f_toread = batch_file_name+str(i)+".csv"
    subprocess.run("gsutil cp "+bucket+"/tmirmira/data/qc/"+f_toread+" ./", shell=True, check=True)
    df = pd.read_csv(f_toread)
    dfs.append(df)

all_ratios = pd.concat(dfs)
all_ratios.to_csv("wgs_ratios_all.csv", index=False)
subprocess.run("gsutil cp wgs_ratios_all.csv "+bucket+"/tmirmira/data/qc", shell=True, check=True)
