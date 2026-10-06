
"""
Builds a maximum-likelihood phylogenetic tree for each group's trimmed alignment using
IQ-TREE. Outputs IQ-TREE result files per group.
"""


import os
from Bio import SeqIO
import subprocess
import glob
import datetime

def build_trees(align_path, iqtree_path, iqtree_conda_env, iqtree_parameters, log, n_threads, min_align_len):

    # Convert from string to integer
    min_align_len = int(min_align_len)

    # Get all alignment file names
    align_files = glob.glob(os.path.join(align_path, "*_trimmed_aligned.fasta"))
    align_files.sort()

    n = 0 # counter
    max_trees = len(align_files) # Max number of trees
    for align in align_files:

        base = os.path.basename(align)
        output_files_pre = os.path.join(iqtree_path, base.replace("_trimmed_aligned.fasta",""))
        
        # Check if alignment is shorter than threshold
        align_len = len(next(SeqIO.parse(align, "fasta")).seq)
        if align_len < min_align_len:
            log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {base} excluded, has {align_len} amino acids in length (<{min_align_len})")
            max_trees -= 1
            continue

        # Run IQ-TREE        
        subprocess.run([iqtree_conda_env,"-s",align,*iqtree_parameters.split(),"-nt","AUTO","-ntmax",n_threads,"--seqtype","AA",
                        "-pre",output_files_pre], check=True)

        n += 1
        if n % 10 == 0:
            log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {n}/{max_trees} trees built.")

    log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : IQ-TREE finished. {n} trees built.")
