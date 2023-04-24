#!/bin/bash
: '
run bash parallel.sh {target_path} {metric}

Summarize summary for given metric and target.

Parameters:
target_path (str): Relative path to target directory
method (str): Desired metric to calculate. Choose from "all", "inf", "clashscore", "gdt", "lddt" or "tm_score"

'

target=$1
metric=$2

# Remove trailing slash from target path if it exists
target=${1%/}

echo "sbatch parallel/parallel.slurm ${target} $metric"
sbatch parallel/parallel_summary.slurm ${target} "$metric"
'