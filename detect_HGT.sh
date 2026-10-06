#!/bin/bash -l

#OAR -l walltime=500:00:00

# Get config file path
config_file="${PWD}/$1"

# Change to script directory
cd "$(dirname "$0")/bin" || exit 1

source /mnt/yeast-data/yglaluno1/miniconda3/etc/profile.d/conda.sh
conda activate python_joaoc

python S00_run_pipeline.py $config_file
