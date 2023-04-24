from target import consolidate_target
import sys

target = sys.argv[1]
metric = sys.argv[2]

# Check if metric is valid
if metric not in ["all", "gdt", "clashscore", "clash", "lddt", "tm_score", "inf"]:
    print("Invalid metric")
    exit(1)

consolidate_target(target, metric=metric)