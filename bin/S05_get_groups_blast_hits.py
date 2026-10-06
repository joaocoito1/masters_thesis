
"""
For each HGT candidate group, builds a FASTA file containing the group's sequences
and their top BLAST hits. Additionally, adds the group's top ingroup (yeast) hits
if the group doesn't already have enough ingroup hits among its top BLAST hits.
"""


import pandas as pd
import os
from Bio import SeqIO
import glob
import datetime

def get_groups_blast_hits(blastp_path, hgt_groups_path, groups_path, AlienIndex_fasta, log, n_top, n_top_ingroup):

    # Convert from string to integer
    n_top = int(n_top)
    n_top_ingroup = int(n_top_ingroup)

    ### Extract BLAST results ###

    # BLAST columns (outfmt 6 format)
    columns = ["query", "subject", "percent_identity", "alignment_length", "mismatches", "gap_opens",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]

    # Extract BLAST file paths
    blast_files = glob.glob(os.path.join(blastp_path, "*blastp.out"))

    # Get all BLAST results
    blast_results = pd.concat([pd.read_csv(f, sep="\t", header=None, names=columns) 
                               for f in blast_files], ignore_index=True)

    # Sort results by query and bitscore
    sorted_blast_results = (blast_results.sort_values(by=["query", "bitscore"],ascending=[True, False])
                                         .drop_duplicates(subset=["query", "subject"], keep="first"))

    # Get top BLAST results
    top_blast_results = (sorted_blast_results
                        .groupby("query")
                        .head(n_top)
                        .reset_index(drop=True))

    # Get top BLAST hits per query -> {query : set of hits}
    query_to_hits = top_blast_results.groupby("query")["subject"].apply(set).to_dict()


    ### Extract groups and BLAST hits ###

    # Get groups
    with open(hgt_groups_path) as f:
        hgt_groups = {}
        for line in f.readlines():
            group_number, group = line.strip().split(": ")
            group = set(group.split("\t"))
            hgt_groups[group_number] = group
        
    # Add top BLAST hits to groups
    hgt_groups_hits = {g_number: set(names) for g_number, names in hgt_groups.items()} # Duplicate dictionary
    for n, group in hgt_groups.items():
        for p in group:
            hgt_groups_hits[n].update(query_to_hits[p])

    ### Add ingroup hits ###

    # Map each query to its group number
    query_to_group = {q: n for n, group in hgt_groups.items() for q in group}

    # Get only ingroup hits
    ingroup_only = (sorted_blast_results[sorted_blast_results["subject"].str.startswith("ingroup")]
                 .assign(group=lambda df: df["query"].map(query_to_group))
                 .dropna(subset=["group"])
                 .sort_values(["group", "bitscore"], ascending=[True, False]))

    # Get top ingroup BLAST hits per group
    top_ingroup_per_group = (ingroup_only
                            .drop_duplicates(subset=["group", "subject"], keep="first")
                            .groupby("group").head(n_top_ingroup)
                            .groupby("group")["subject"].apply(list).to_dict())

    # Add ingroup hits to groups that need it
    for n, group in hgt_groups.items():
        ingroup_count = sum(1 for seq_id in hgt_groups_hits[n] if seq_id.startswith("ingroup"))
        
        if ingroup_count < n_top_ingroup:
            hgt_groups_hits[n].update(top_ingroup_per_group.get(n, []))

    # Get all sequences used
    fasta_blast_db = SeqIO.index(AlienIndex_fasta, "fasta")

    # Create a FASTA file with each group
    for n, group in hgt_groups_hits.items():

        records = [fasta_blast_db[seq_id] for seq_id in group]

        fasta_file = os.path.join(groups_path, f"{n}.fasta")
        SeqIO.write(records, fasta_file, "fasta")


    log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : Top BLAST hits added to groups.")
