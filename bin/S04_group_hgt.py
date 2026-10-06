
"""
Identifies all HGT candidates (Alien Index above threshold) across strains and groups
them based on direct BLAST hits between candidates (two candidates are placed in the
same group if one is a BLAST hit of the other, with groups merging if they share members).
Outputs a single file listing candidate groups (one group per line).
"""


import pandas as pd
import os
import datetime

def group_hgt(strains_df, AI_path, blastp_path, hgt_groups_path, log, AlienIndex_threshold):

    # Extract strain names to be analyzed
    strain_names = strains_df["strain_name"].tolist()

    ##### Extract candidate HGT sequence names #####
    candidate_hgt = {}
    for strain in strain_names:

        # Build input file path
        AI_file = os.path.join(AI_path, f"{strain}_AlienIndex_results.csv")

        AlienIndex_results = pd.read_csv(AI_file)

        # Extract all sequences with Alien Index above threshold
        candidate_hgt[strain] = set(AlienIndex_results.query(f"AI>{AlienIndex_threshold}")["query"])

    # Get all HGT candidates in one set
    all_candidate_hgt = set().union(*candidate_hgt.values())


    ##### Group candidate HGT #####

    # BLAST columns (outfmt 6 format)
    columns = ["query", "subject", "percent_identity", "alignment_length", "mismatches", "gap_opens",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]

    # Extract BLAST file paths
    blast_files = [os.path.join(blastp_path, f"{s}_AlienIndex_blastp.out") 
                   for s in strain_names]

    # Get all BLAST results (if files exist)
    blast_results = pd.concat([pd.read_csv(f, sep="\t", header=None, names=columns) 
                               for f in blast_files if os.path.exists(f)], ignore_index=True)

    # Get BLAST hits per query -> {query : set of hits}
    query_to_hits = blast_results.groupby("query")["subject"].apply(set).to_dict()

    hgt_groups = [] # will contain sets, each corresponding to a group
    grouped_proteins = set() # will contain already grouped HGT candidates
    for prot in all_candidate_hgt:

        # Get BLAST hits that are other HGT candidates
        prot_hit_hgt = {p for p in query_to_hits[prot] if p in all_candidate_hgt}
        prot_hit_hgt.add(prot)
        
        # Get ones already in a group
        in_group = {p for p in prot_hit_hgt if p in grouped_proteins}

        # Add to the set the ones that will be grouped bellow
        grouped_proteins.update(prot_hit_hgt)

        # Get group numbers
        in_group_numbers = {g_numb 
                            for g_numb, elems in enumerate(hgt_groups) 
                            for p in in_group 
                            if p in elems}
        
        # Number of different groups
        numb_groups = len(in_group_numbers)

        # If none already in group, create new one
        if numb_groups == 0:
            hgt_groups.append(prot_hit_hgt)
        
        else:
            # If any already in a group, add to that
            if numb_groups == 1:
                group_number = in_group_numbers.pop()
                hgt_groups[group_number].update(prot_hit_hgt)

            # If already in multiple groups, join all groups and add to that group
            elif numb_groups > 1:
                joint_group_number = min(in_group_numbers)

                # Add all elements of other intersecting groups to joint group
                in_group_numbers.remove(joint_group_number)
                for n in in_group_numbers:
                    # Add elements to intersecting group
                    hgt_groups[joint_group_number].update(hgt_groups[n])
                    # Mark group as deleted (avoid index shifting)
                    hgt_groups[n] = None
                
                hgt_groups[joint_group_number].update(prot_hit_hgt)
        
        # Remove deleted groups
        hgt_groups = [g for g in hgt_groups if g is not None]

    # Sort by group size (descending)
    hgt_groups.sort(key=len, reverse=True)

    # Save groups to file
    with open(hgt_groups_path, "w") as f:
        for i, hgt_candidates in enumerate(hgt_groups):
            f.write(f"G{i:07d}: {'\t'.join(hgt_candidates)}\n")

    log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {len(all_candidate_hgt)} HGT candidates in {len(hgt_groups)} groups.")
    log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {os.path.basename(hgt_groups_path)} created.")

    return all_candidate_hgt