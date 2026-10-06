#!/bin/bash
set -e -o pipefail

# Runs BLASTP for each strain's proteome against the AlienIndex DB
# (bacteria/fungi/yeast/WS clade) to gather hit data used for Alien
# Index calculation in the next step. One BLAST output file per strain.


##############
### Inputs ###
##############

# Path to run blastp
blastp_path="$1"

# Path to csv with strain information
strains_csv="$2"

# BLAST Database path + name
AlienIndex_BLAST_DB="$3"

# Path to BLAST output folder
cd "$4"

# Path to log file
log_file="$5"

# Number of threads
n_threads="$6"

# BLAST e-value threshold
evalue="$7"

# BLAST maximum number of target sequences
max_seqs="$8"


##############
###  Main  ###
##############


while IFS=',' read -r strain proteome; do
    
    $blastp_path -db $AlienIndex_BLAST_DB -query $proteome -out ${strain}_AlienIndex_blastp.out \
        -outfmt 6 -num_threads $n_threads -evalue $evalue -max_target_seqs $max_seqs

    message="$(date "+%Y-%m-%d %H:%M:%S") : ${strain}_AlienIndex_blastp.out created."
    echo "$message"
    echo "$message" >> "$log_file"

done < <(tail -n +2 "$strains_csv") # skip header