#!/bin/bash
set -e -o pipefail

# Aligns each HGT candidate group's FASTA file with MAFFT, skipping groups below
# a minimum sequence count, then trims each resulting alignment with trimAL.
# Outputs one aligned and trimmed FASTA file per group.


##############
### Inputs ###
##############

# Path to run MAFFT
mafft_path="$1"

# Path to run trimAL
trimal_path="$2"

# Path to folder with FASTA groups
cd "$3"

# Path to aligments output folder
align_path="$4"

# Path to log file
log_file="$5"

# Number of threads
n_threads="$6"

# Minimum number of sequences per group
min_numb_seqs="$7"

# MAFFT algorith parameters
mafft_parameters="$8"

# trimAL trimming method
trimal_method="$9"


##############
###  Main  ###
##############

mafft_count=0 #count
total_groups=$(ls -1 | wc -l)

### MAFFT ###
for group in *; do

    # Number of sequences in group
    num_seqs=$(grep -c '^>' $group)

    if (( num_seqs < min_numb_seqs )); then
        message="$(date "+%Y-%m-%d %H:%M:%S") : $group does not have enough sequences (it has $num_seqs)"
        echo "$message"
        echo "$message" >> "$log_file"

    else
        # Output file path (G0000000_aligned.fasta)
        mafft_file="${align_path}/${group%.fasta}_aligned.fasta"

        # Run MAFFT
        $mafft_path $mafft_parameters --thread $n_threads $group > $mafft_file

        # Print message every 25 alignments
        mafft_count=$((mafft_count + 1))
        if (( mafft_count % 25 == 0 )); then
            message="$(date "+%Y-%m-%d %H:%M:%S") : $mafft_count/$total_groups groups aligned."
            echo "$message"
            echo "$message" >> "$log_file"
        fi
    fi
done

message=$'\n'"$(date "+%Y-%m-%d %H:%M:%S") : MAFFT finished. $mafft_count groups aligned."$'\n'
echo "$message"
echo "$message" >> "$log_file"


# Move to alignments folder
cd $align_path


### trimAL ###
for mafft_file in *_aligned.fasta; do

    # trimAL output file path
    trimal_file="${align_path}/${mafft_file%_aligned.fasta}_trimmed_aligned.fasta"    

    # Run trimAL
    $trimal_path -in "$mafft_file" -out "$trimal_file" -"${trimal_method}"

done

message="$(date "+%Y-%m-%d %H:%M:%S") : trimAL finished."$'\n'
echo "$message"
echo "$message" >> "$log_file"
