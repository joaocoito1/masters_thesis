
"""
Calculates the Alien Index (AI) for each candidate protein per strain, based on
normalized best bitscore comparison between outgroup (bacteria/fungi) and ingroup
(non-W/S yeast) BLAST hits. Outputs one CSV per strain with AI scores and top hit subjects.
"""


import pandas as pd
import os
import datetime

def calculate_AlienIndex(strains_df, blastp_path, AI_path, log, AlienIndex_threshold):

    # Extract strain names to calculate Alien Index
    strain_names = strains_df["strain_name"].tolist()

    # BLAST columns (outfmt 6 format)
    columns = ["query", "subject", "percent_identity", "alignment_length", "mismatches", "gap_opens",
            "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]

    for strain in strain_names:

        # Build file path
        strain_blast_out = os.path.join(blastp_path, f"{strain}_AlienIndex_blastp.out")

        # Load file to dataframe
        blast_df = pd.read_csv(strain_blast_out, sep="\t", header=None, names=columns)

        # Sort results by query and bitscore
        blast_df = blast_df.sort_values(by=["query", "bitscore"],ascending=[True, False])

        # Get all queries
        all_queries = sorted(set(blast_df["query"]))


        # Compute best bitscores (no hit -> 0)
        outgroup_scores = (blast_df[blast_df["subject"].str.startswith(("outgroup"))]
                            .groupby("query")["bitscore"].max()
                            .reindex(all_queries, fill_value=0))

        ingroup_scores = (blast_df[blast_df["subject"].str.startswith("ingroup")]
                            .groupby("query")["bitscore"].max()
                            .reindex(all_queries, fill_value=0))

        self_scores = (blast_df[blast_df["query"] == blast_df["subject"]]
                        .groupby("query")["bitscore"].max()
                        .reindex(all_queries, fill_value=0))
        
        # Normalization
        outgroup_norm_scores = outgroup_scores / self_scores
        ingroup_norm_scores = ingroup_scores / self_scores


        # Get subjects
        outgroup_subject = (blast_df[blast_df["subject"].str.startswith(("outgroup"))]
                            .groupby("query")["subject"].first()
                            .reindex(all_queries, fill_value="No Subject"))
        
        ingroup_subject  = (blast_df[blast_df["subject"].str.startswith("ingroup")]
                            .groupby("query")["subject"].first()
                            .reindex(all_queries, fill_value="No Subject"))

        # Create DataFrame
        results = pd.DataFrame({
            "outgroup_norm_bitscore": outgroup_norm_scores,
            "ingroup_norm_bitscore": ingroup_norm_scores,
            "self_bitscore": self_scores,
            "outgroup_subject": outgroup_subject,
            "ingroup_subject": ingroup_subject,
        })

        # Turn the index (query IDs) into a real column
        results = results.reset_index().rename(columns={"index": "query"})

        # Calculate Alien Index
        results["AI"] = results["outgroup_norm_bitscore"] - results["ingroup_norm_bitscore"]

        # Get number of HGT candidates
        numb_candidate_hgt = len(results.query(f" AI > {AlienIndex_threshold}"))

        # Save output
        output_file = os.path.join(AI_path, f"{strain}_AlienIndex_results.csv")
        results.to_csv(output_file, index=False)
        log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {strain}_AlienIndex_results.csv created. {numb_candidate_hgt} HGT candidates (AI > {AlienIndex_threshold})")