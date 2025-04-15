#!/bin/bash
#SBATCH -A standby                 # or your account/queue
#SBATCH --nodes=1                  # 1 node
#SBATCH --ntasks=100               # use 100 CPU cores on that node
#SBATCH --time=4:00:00            # e.g., 72 hours max
#SBATCH --job-name=dapper_run
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

# 1) Parse arguments from the command line
iter=$1
n=$2
d=$3
p=$4
mu=$5
sa=$6

# 2) Create an output directory named by parameter settings
outdir="dapper_comparison_output_${iter}_${n}_${d}_${p}_${mu}_${sa}"
mkdir -p "$outdir"

# 3) Load modules (or activate environment) needed for R
module load R/4.4.1

# 4) Move to the submission directory
cd "$SLURM_SUBMIT_DIR"

echo "Running Dapper comparison with iter=$iter n=$n d=$d p=$p mu=$mu sa=$sa"
echo "Output folder: $outdir"

# 5) Run your R script
Rscript dapper_comparison.R $iter $n $d $p $mu $sa

# 6) Move results and SLURM logs into the output folder
mv sim_results.RData "$outdir/" 2>/dev/null
mv slurm_${SLURM_JOB_ID}.out "$outdir/" 2>/dev/null
mv slurm_${SLURM_JOB_ID}.err "$outdir/" 2>/dev/null
