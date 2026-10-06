

###################
###   Imports   ###
###################

from pipeline_utils import *

from S02_calculate_AlienIndex import calculate_AlienIndex
from S03_blast_missing_hits import blast_missing_hits
from S04_group_hgt import group_hgt
from S05_get_groups_blast_hits import get_groups_blast_hits
from S07_build_trees import build_trees
from S09_generate_review_excel import generate_review_excel

import subprocess
import os
import pandas as pd
import sys
import datetime


####################
###    Inputs    ###
####################

## Get provived config file path ##
config_file_path = sys.argv[1]

## Extract inputs from configs file ##
configs = extract_inputs(config_file_path)

## Define the number of threads ##
n_threads = configs["Parameters"]["n_threads"] if "n_threads" in configs["Parameters"] else str(os.cpu_count())

## Check if all the necessary inputs are provided ##
check_input(configs)
print("Input check succeeded!")

## Create Run folder ##
run_path, path = create_run_folder(configs, config_file_path)

## Extract strain information ##
strains_df = pd.read_csv(configs["Paths"]["strains_csv"])
strains_df = strains_df.set_index("strain_name",drop=False)

## Print initial informations ##
log = make_logger(path["log_file"])
log(print_message["00"])
initial_prints(log, run_path, path, strains_df, n_threads, configs)

####################
###   Pipeline   ###
####################

## 01 - Alien Index BLAST ##
log(print_message["01"])
os.makedirs(path["blast"])
subprocess.run(["bash", "S01_blast_AlienIndex.sh", configs["Tool Paths"]["blastp_path"], configs["Paths"]["strains_csv"], 
                configs["Paths"]["AlienIndex_BLAST_DB"], path['blast'], path['log_file'], n_threads, 
                configs["Parameters"]["evalue"], configs["Parameters"]["max_seqs"]], check=True)


## 02 - Calculate Alien Index ##
log(print_message["02"])
os.makedirs(path["AI"])
calculate_AlienIndex(strains_df, path['blast'], path['AI'], log, configs["Parameters"]["AlienIndex_threshold"])


## 03 - BLAST HGT candidates vs Ingroup DB ##
log(print_message["03"])
blast_missing_hits(strains_df, path['AI'], path['blast'], configs["Tool Paths"]["blastp_path"], configs["Paths"]["Ingroup_BLAST_DB"], 
                   log, n_threads, configs["Parameters"]["AlienIndex_threshold"], configs["Parameters"]["evalue_2nd_blast"], 
                   configs["Parameters"]["max_seqs_2nd_blast"])


## 04 - Group HGT candidates ##
log(print_message["04"])
all_candidate_hgt = group_hgt(strains_df, path['AI'], path['blast'], path["hgt_groups"], log,
                              configs["Parameters"]["AlienIndex_threshold"])


## 05 - Add top BLAST hits to groups ##
log(print_message["05"])
os.makedirs(path["groups"])
get_groups_blast_hits(path['blast'], path['hgt_groups'], path['groups'], configs["Paths"]["AlienIndex_fasta"], log,
                      configs["Parameters"]["n_top_hits"], configs["Parameters"]["n_top_ingroup_hits"])


## 06 - Align sequences by group ##
log(print_message["06"])
os.makedirs(path["align"])
subprocess.run(["bash", "S06_align_groups.sh", configs["Tool Paths"]["mafft_path"], configs["Tool Paths"]["trimal_path"], 
                path['groups'], path['align'], path['log_file'], n_threads, 
                configs["Parameters"]["min_numb_seqs"], configs["Parameters"]["mafft_parameters"], configs["Parameters"]["trimal_method"]], 
                check=True)


## 07 - Build a phylogenetic tree for each group ##
log(print_message["07"])
os.makedirs(path["iqtree"])
build_trees(path['align'], path['iqtree'], configs["Tool Paths"]["iqtree_path"], configs["Parameters"]["iqtree_parameters"], 
            log, n_threads, configs["Parameters"]["min_align_len"])


## 08 - Root and color trees ##
log(print_message["08"])
os.makedirs(path["trees"])
subprocess.run([configs["Tool Paths"]["R_path"], "S08_edit_tree.R", path['groups'], path['iqtree'], path['trees'], path['groups_wo_trees'],
                configs["Parameters"]["names_groups"], configs["Parameters"]["color_groups"], configs["Parameters"]["recipient_color"],
                path['hgt_groups'], path['log_file']], check=True)


## 09 - Generate excel to manually classify HGT candidates ##
log(print_message["09"])
generate_review_excel(strains_df, log, path['AI'], path['hgt_groups'], path['groups_wo_trees'], path['hgt_review'])

## Final message ##
log(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
log(print_message["end"])
