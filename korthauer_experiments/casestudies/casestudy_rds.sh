#!/bin/bash
#SBATCH --cpus-per-task 1
R CMD BATCH --no-save casestudies.R
