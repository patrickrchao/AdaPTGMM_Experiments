#!/bin/bash
#SBATCH --cpus-per-task 14
R CMD BATCH --no-save run_all_sims.R
