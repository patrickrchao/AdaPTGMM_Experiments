#!/bin/bash
#SBATCH --cpus-per-task 13
R CMD BATCH --no-save generate_alpha_summary.R
