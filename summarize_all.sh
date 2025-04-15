#!/bin/bash
# summarize_all.sh
#
# This script automatically finds all simulation output files
# with the directory prefix "/home/navarr72/dpwork/dapper_comparison_output_"
# and then calls summarize_results.R with those files as arguments.
#
# The expected file structure is:
#   /home/navarr72/dpwork/dapper_comparison_output_<params>/sim_results.RData

# Use a bash glob to find all sim_results.RData files
files=(/home/navarr72/dpwork/dapper_comparison_output_*/sim_results.RData)

# Optionally, you can echo the files to verify
echo "Found ${#files[@]} simulation output files:"
for f in "${files[@]}"; do
    echo "$f"
done

# Call the R script with all file names as arguments
Rscript summarize_results.R "${files[@]}"
