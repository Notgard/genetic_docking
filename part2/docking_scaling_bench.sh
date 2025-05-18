#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=03:00:00
#SBATCH --mem=20G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name="docking_scaling_test"
#SBATCH --error=err/scaling.err
#SBATCH --output=out/scaling.out

romeo_load_x64cpu_env

# Clean previous results
rm -f results/*
rm -f fitness_log.csv

# Paths
LIGAND="../files/ligand_NO_BOND.mol2"
SITE="../files/site.mol2"
EXEC="./genetic_algorithm"

# Header
echo "THREADS,CPU_TIME,OMP_WALLTIME" > scaling_results.csv

# Function to extract timing results from program output
extract_timing() {
    CPU=$(grep "Fortran CPU time elapsed" $1 | awk '{print $5}')
    OMP=$(grep "OpenMP Walltime elapsed" $1 | awk '{print $5}')
    echo "$2,$CPU,$OMP" >> scaling_results.csv
}

# --- 1. Sequential test
export OMP_NUM_THREADS=1
echo "Running sequential version (1 thread)..."
OUTFILE="out/sequential.out"
$EXEC $LIGAND $SITE > $OUTFILE
extract_timing $OUTFILE 1

# --- 2. Parallel scaling tests
for t in 2 4 8 16 32 64; do
    export OMP_NUM_THREADS=$t
    echo "Running OpenMP version with $t threads..."
    OUTFILE="out/omp_${t}threads.out"
    $EXEC $LIGAND $SITE > $OUTFILE
    extract_timing $OUTFILE $t
done

echo "Scaling test complete. Results saved to scaling_results.csv"
