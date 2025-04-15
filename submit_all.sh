#!/bin/bash

# Submit all combos in a single pass. Each sbatch call submits a separate job
# that runs run_dapper.sh with the arguments we pass in.

# 1) n-values
for n in 100 200 400 1600 3200; do
  sbatch run_dapper.sh 5000 $n 4 1 4 2
done

# 2) d-values
for d in 2 8 16 32 64; do
  sbatch run_dapper.sh 5000 800 $d 1 4 2
done

# 3) p-values
for p in 0.01 0.05 0.1; do
  sbatch run_dapper.sh 5000 800 4 $p 4 2
done

# 4) mu-values
for mu in 1 2 8; do
  sbatch run_dapper.sh 5000 800 4 1 $mu 2
done

# 5) sa-values
for sa in 1 4 8; do
  sbatch run_dapper.sh 5000 800 4 1 4 $sa
done

# 6) Finally run the default combo once
sbatch run_dapper.sh 5000 800 4 1 4 2
