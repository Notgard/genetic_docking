#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=03:00:00
#SBATCH --mem=20G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name "genetic_docking"
#SBATCH --error=err/job.err
#SBATCH --output=out/job.out

romeo_load_x64cpu_env

./genetic_algorithm ../files/ligand_NO_BOND.mol2 ../files/site.mol2
