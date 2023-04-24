#!/bin/bash
: '
run bash parallel.sh {target_path} {metric}

Computing cluster script to find all metrics for a given target.

Parameters:
target_path (str): Relative path to target directory
method (str): Desired metric to calculate. Choose from "all", "inf", "clashscore", "gdt", "lddt" or "tm_score"
'


target=$1
metric=$2

# echo "r1107 r1108 r1116 r1117 r1126 r1128 r1136 r1138 r1149 r1156 r1189 r1190" | xargs -n 1 -P 10 sh -c '

# Remove trailing slash from target path if it exists
target=${1%/}

# Check if metric is "all"
if [ "$metric" = "all" ]; then
    # Run all metrics
    echo "sbatch parallel/parallel.slurm ${target} inf"
    sbatch parallel/parallel.slurm ${target} "inf"
    
    echo "sbatch parallel/parallel.slurm ${target} clashscore"
    sbatch parallel/parallel.slurm ${target} "clashscore"

    echo "sbatch parallel/parallel.slurm ${target} gdt"
    sbatch parallel/parallel.slurm ${target} "gdt"

    echo "sbatch parallel/parallel.slurm ${target} lddt"
    sbatch parallel/parallel.slurm ${target} "lddt"

    echo "sbatch parallel/parallel.slurm ${target} tm_score"
    sbatch parallel/parallel.slurm ${target} "tm_score"

else
    # Run only the specified metric
    echo "sbatch parallel/parallel.slurm ${target} $metric"
    sbatch parallel/parallel.slurm ${target} "$metric"
fi
