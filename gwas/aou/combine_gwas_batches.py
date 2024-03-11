import pandas as pd
import subprocess
import os
import sys

bucket = os.getenv('WORKSPACE_BUCKET')

phenotype = sys.argv[1]
gwas_file_prefix = sys.argv[2] #e.g. ldl_cholesterol_hail_EUR_WHITE

chr_regions_file = "chrom_regions_hg38.txt"
with open(chr_regions_file, "r") as f:
    chr_regions = [line.strip() for line in f.readlines()]

dfs = []
total_snps = 0
comment_lines = None
for r in chr_regions:
    chrom = r.split(":")[0]
    full_file = gwas_file_prefix + "_%s"%(r.replace(":", "_").replace("-","_")) + ".gwas.tab"
    subprocess.run("gsutil cp "+bucket+"/gwas/"+phenotype+"/"+full_file+" ./", shell=True, check=True)
    if chrom == "chr1":
        #read header lines
        with open(full_file, 'r') as f:
            comment_lines = [line.strip() for line in f if line.startswith('#')]
    df = pd.read_csv(full_file, sep="\t", comment="#")
    total_snps += df.shape[0]
    dfs.append(df)


first_line = comment_lines[0].split(" ")
new_first_line = []
skip = []
for i in range(len(first_line)):
    if i in skip:
        continue
    if first_line[i] == "--region":
        skip.append(i+1)
    else:
        new_first_line.append(first_line[i])
new_first_line = " ".join(new_first_line)
comment_lines[0] = new_first_line


full_gwas = pd.concat(dfs)
assert(full_gwas.shape[0] == total_snps)
full_gwas.columns = dfs[0].columns
full_gwas_file = gwas_file_prefix + ".gwas.tab"
with open(full_gwas_file, 'w') as f:
    for comment in comment_lines:
        f.write(comment + '\n')
    full_gwas.to_csv(f, sep="\t", index=False, header=True)

subprocess.run("gsutil cp "+full_gwas_file+" "+bucket+"/gwas/"+phenotype+"/"+full_gwas_file, shell=True, check=True)
