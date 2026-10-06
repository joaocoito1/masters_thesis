
import os
import datetime
import shutil
import pandas as pd
import datetime
import sys

########################################
###  File/Folder location and names  ###
########################################

dir_name = {"base":"",
            "intermediate_results" : "Intermediate_Results",
            "final_results" : "Final_Results",
            "blast" : "BLAST",
            "groups" : "Groups",
            "align" : "Alignments",
            "iqtree" : "IQ-TREE",
            "AI" : "AlienIndex",
            "trees" : "Trees"}

dir_location = {"blast" : "intermediate_results",
                "groups" : "intermediate_results",
                "align" : "intermediate_results",
                "iqtree" : "intermediate_results",
                "AI" : "final_results",
                "trees" : "final_results"}


file_name = {"log_file":"Run.log",
             "config_file":"config.ini",
             "strains_csv":"strains.csv",
             "hgt_groups":"Groups_candidates_HGT.txt",
             "groups_wo_trees":"Groups_without_trees.txt",
             "hgt_review":"HGT_review.xlsx"}

file_location = {"log_file":"base",
                 "config_file":"base",
                 "strains_csv":"base",
                 "hgt_groups":"final_results",
                 "groups_wo_trees":"trees",
                 "hgt_review":"final_results"}



relative_path = {folder_name : os.path.join(dir_name[dir_location[folder_name]], dir_name[folder_name]) 
                 for folder_name in dir_location}

relative_path.update({f_name : os.path.join(dir_name[file_location[f_name]], file_name[f_name])
                               if file_location[f_name] not in dir_location else
                               os.path.join(relative_path[file_location[f_name]], file_name[f_name])
                      for f_name in file_location})


#######################################
###       Auxiliary functions       ###
#######################################

def extract_inputs(config_file):
    """ 
    Extracts information from a .ini file into a dictionary of dictionaries.
    Comments (Text after "#") are ignored
    """

    with open(config_file) as f:
        lines = f.readlines()

    configs = {}
    for l in lines:
        l = l.split("#")[0].strip() # remove everything after "#" (comments)

        if "[" in l:
            configs[l[1:-1]] = {}
            current_config = l[1:-1]
        
        elif "=" in l:
            name, value = l.split("=")
            configs[current_config][name.strip()] = value.strip()
    
    return configs



def create_run_folder(configs, config_file_path, relative_path=relative_path):

    ## Define Run output path and create it ##
    if "run_name" in configs["Names"]:
        run_name = configs["Names"]["run_name"]
    else:
        run_name = datetime.datetime.now().strftime("%Y%b%d")
        
    run_path = os.path.join(configs["Paths"]["output_path"], "Results_"+run_name)
    
    while os.path.exists(run_path):
        run_path += "_1"
    run_path = os.path.abspath(run_path)

    ## Define folder/file absolute paths ##
    absolute_paths = {name : os.path.join(run_path, relative_path[name]) 
            for name in relative_path}
    
    ## Create file directories ##
    for f in ["log_file","config_file","strains_csv"]:
        os.makedirs(os.path.dirname(absolute_paths[f]), exist_ok=True)

    ## Copy config and strains files ##
    shutil.copy(config_file_path, absolute_paths["config_file"])
    shutil.copy(configs["Paths"]["strains_csv"], absolute_paths["strains_csv"])
    
    ## Create run file ##
    open(absolute_paths["log_file"], "w").close()

    
    return run_path, absolute_paths # absolute_paths = {name : path}


#######################################
###          Print messages         ###
#######################################

# Nested function
def make_logger(log_file):
    """ Returns a function that prints a message and appends it to log_file """
    def log(*args, sep=" "):
        """ Prints a message and saves it to log file """
        message = sep.join(str(a) for a in args)
        print(message)
        with open(log_file, "a") as lf:
            lf.write(message + "\n")
    return log



def initial_prints(log, run_path, path, strains_df, n_threads, configs, size=25):

    log(f"{'Run:':<{size}}{run_path}")
    log(f"{'Log file:':<{size}}{path['log_file']}")
    log(f"{'Config file:':<{size}}{path['config_file']}")
    log(f"{'Started:':<{size}}{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    strains_label = f"Strains ({len(strains_df['strain_name'])}):"
    log(f"{strains_label:<{size}}{', '.join(strains_df['strain_name'])}")

    log(f"{'Number of threads:':<{size}}{n_threads}")
    log(f"{'Alien Index threshold:':<{size}}{configs['Parameters']['AlienIndex_threshold']}")



print_message = {"00":"""
###############################################################
##                                                           ##
##         HGT Detection Pipeline for the W/S Clade          ##
##                                                           ##
###############################################################
""",
"01":"""\n
###############################################################
##                  01 - Alien Index BLAST                   ##
###############################################################
""",
"02":"""\n
###############################################################
##                 02 - Calculate Alien Index                ##
###############################################################
""",
"03":"""\n
###############################################################
##          03 - BLAST HGT candidates vs Ingroup DB          ##
###############################################################
""",
"04":"""\n
###############################################################
##                 04 - Group HGT candidates                 ##
###############################################################
""",
"05":"""\n
###############################################################
##             05 - Add top BLAST hits to groups             ##
###############################################################
""",
"06":"""\n
###############################################################
##               06 - Align sequences by group               ##
###############################################################
""",
"07":"""\n
###############################################################
##       07 - Build a phylogenetic tree for each group       ##
###############################################################
""",
"08":"""\n
###############################################################
##                 08 - Root and color trees                 ##
###############################################################
""",
"09":"""\n
###############################################################
##  09 - Generate excel to manually classify HGT candidates  ##
###############################################################
""",
"end":"""\n
###############################################################
##                     Pipeline finished!                    ##
###############################################################
"""}




def get_fasta_ids(fasta_path):
    """ Extract sequence IDs from a FASTA file """
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0].strip())  # header up to first space
                if "," in line:
                    sys.exit("There are commas (,) in the headers.")
    return ids



def check_input(configs):
    """
    Check if all the necessary inputs are correct
    """

    ## Missing files and parameters ##

    errors = []

    configs_requirements = {'Paths': {'AlienIndex_fasta', 'AlienIndex_BLAST_DB', 'strains_csv', 'output_path', 'Ingroup_BLAST_DB'}, 
                            'Parameters': {'n_top_hits', 'trimal_method', 'max_seqs', 'evalue', 'color_groups', 'min_numb_seqs', 
                                           'names_groups', 'min_align_len', 'AlienIndex_threshold', 'n_top_ingroup_hits', 'mafft_parameters', 
                                           'recipient_color', 'evalue_2nd_blast', 'iqtree_parameters', 'max_seqs_2nd_blast'}, 
                            'Tool Paths': {'blastp_path', 'trimal_path', 'mafft_path', 'iqtree_path', 'R_path'}}

    for section in configs_requirements:
        if section not in configs:
            errors.append(f"ERROR: Missing section [{section}] from config file.")

        else:
            # If section exists
            for name in configs_requirements[section]:
                if name not in configs[section]:
                    errors.append(f"ERROR: {name} not in section [{section}] from config file.")

                else:
                    # If name is in config
                    x = configs[section][name]

                    if section in ["Tool Paths","Paths"]:
                        if name in ["AlienIndex_BLAST_DB","Ingroup_BLAST_DB"]:
                            x += ".pdb" # if is BLAST db, add extension to search for file
                        if not os.path.exists(x):
                            errors.append(f"ERROR: {name} file not found: {x}")
    # Print errors
    if errors:
        sys.exit("Input check failed:\n- " + "\n- ".join(errors))

    print("Necessary files are present.")


    ## Check Alien Index database ##

    AI_ids = get_fasta_ids(configs["Paths"]["AlienIndex_fasta"])

    print("Alien Index FASTA sequence IDs imported.")
    
    missing = {}

    strains_df = pd.read_csv(configs["Paths"]["strains_csv"])
    if "strain_name" not in list(strains_df) or "proteome_path" not in list(strains_df):
        sys.exit(f"Missing strain_name or proteome_path: {configs['Paths']['strains_csv']}")


    for strain, proteome_path in zip(strains_df["strain_name"], strains_df["proteome_path"]):
        proteome_ids = get_fasta_ids(proteome_path)
        not_found = proteome_ids - AI_ids
        if not_found:
            missing[strain] = not_found

    if len(missing) > 0:
        sys.exit("Missing proteins in Alien Index FASTA file for the following strains:\n- " + '\n- '.join(missing.keys()) + 
                 "\n\n\n" + '\n\n'.join([f"Strain: {s}\n{m}" for s, m in missing.items()]))

        

def check_proteomes_in_alienindex(strains_df, alienindex_fasta_path):
    alien_ids = get_fasta_ids(alienindex_fasta_path)

    missing = {}
    for strain, proteome_path in zip(strains_df["strain_name"], strains_df["proteome_path"]):
        proteome_ids = get_fasta_ids(proteome_path)
        not_found = proteome_ids - alien_ids
        if not_found:
            missing[strain] = not_found

    return missing
