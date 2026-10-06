
"""
For candidates with AI above threshold and no ingroup (yeast) hit, BLAST their
sequences against a yeast (Saccharomycotina) DB to get hits missed by the first 
BLAST (AlienIndex DB). Outputs one FASTA + BLAST result file per strain.
"""


import pandas as pd
import os
from Bio import SeqIO
import subprocess
import datetime



def blast_missing_hits(strains_df,AI_path,blastp_path,blast_conda_env,Ingroup_BLAST_DB,log,n_threads,
                       AlienIndex_threshold,evalue_2nd_blast,max_seqs_2nd_blast):

    # Extract strain names to be analyzed
    strain_names = strains_df["strain_name"].tolist()

    # Extract candidate HGT sequences without ingroup hits
    hgt_no_ingroup_hits = {}
    for strain in strain_names:
        # Build input file path
        AI_file = os.path.join(AI_path, f"{strain}_AlienIndex_results.csv")

        AlienIndex_results = pd.read_csv(AI_file)

        # Extract all sequences with Alien Index above threshold and no ingroup BLAST hits
        hgt_no_ingroup_hits[strain] = set(AlienIndex_results
                                         .query(f" AI > {AlienIndex_threshold} & ingroup_subject == 'No Subject' ")["query"])


    # Create a FASTA file for each strain with the seqs missing ingroup hits
    for strain, prots in hgt_no_ingroup_hits.items():

        # Number of proteins missing Ingroup
        numb_missing_prots = len(prots)
        
        # If there are no proteins missing Ingroup
        if numb_missing_prots == 0:
            log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {strain} has no candidates missing ingroup hit, skipping BLAST.")
            continue

        log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {strain} has {numb_missing_prots} candidates missing ingroup hit.")

        # Create a FASTA file for each strain with the seqs missing ingroup hits
        
        proteome_path = strains_df.loc[strain, "proteome_path"]

        idx_fasta_proteome = SeqIO.index(proteome_path, "fasta")
        
        records = [idx_fasta_proteome[seq_id] for seq_id in sorted(prots)]

        fasta_file = os.path.join(blastp_path, f"{strain}_HGT_MissingIngroup.fasta")
        SeqIO.write(records, fasta_file, "fasta")


        # Run BLAST

        out_file = os.path.join(blastp_path, f"{strain}_AlienIndex_HGT_MissingIngroup_blastp.out")

        subprocess.run([blast_conda_env,"-db",Ingroup_BLAST_DB,"-query",fasta_file,"-out",out_file,"-outfmt","6",
                        "-num_threads",n_threads,"-evalue",evalue_2nd_blast,"-max_target_seqs",max_seqs_2nd_blast], check=True)

        log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {strain}_AlienIndex_HGT_MissingIngroup_blastp.out created.")