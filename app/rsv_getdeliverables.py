#!/home/sbsuser/miniconda3/bin/python

import os, sys, glob, re, shutil
from pathlib import Path

program_name = Path(__file__).name
project_name = sys.argv[-1]

include_rawreads = False

dir_project = Path("/mnt/tank/bench/projects") / project_name
os.chdir(dir_project)
try:
    shutil.rmtree("Customer_data")
except:
    print("[warning] all existing files from Customer_data will be removed.")
os.mkdir("Customer_data")

def extract_sv(fullname):
    return re.search('.*(SV\d{5}[\-\_]\d{4}).*', fullname).group(1)

def move_files(expr, dest):
    source = glob.glob(expr)
    for file in source:
        shutil.copy2(file, dest/Path(file).name)

print("[rv-g] Getting files for Customer Data Package")
sample_dirs = sorted(glob.glob("RSV_OUTPUT/*"))
sample_names = [extract_sv(d) for d in sample_dirs]
samples = zip(sample_names, sample_dirs)

for spl_name, spl_dir in samples:
    print("Getting files from dir:", spl_dir)
    dir_dest = Path("Customer_data") / spl_name
    os.mkdir(dir_dest)

    # fastq files
    # move_files(f"{spl_dir}/*_trimmed.fastq.gz", dir_dest)
    # move_files(f"{spl_dir}/*_deduped.fastq.gz", dir_dest)
    # move_files(f"{spl_dir}/*_umiextracted.fastq.gz", dir_dest)

    move_files(f"{spl_dir}/*_genome_all.sam", dir_dest) # alignment files
    move_files(f"{spl_dir}/*.vcf", dir_dest) # VCF files
    move_files(f"{spl_dir}/*.txt", dir_dest) # Get summary charts

print("[rv-g] Program complete")
