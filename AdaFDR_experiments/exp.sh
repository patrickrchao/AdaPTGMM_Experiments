#!/bin/bash
#SBATCH --cpus-per-task 1
R CMD BATCH --no-save run_all_exp.R
