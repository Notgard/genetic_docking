#!/bin/bash

# Usage: ./plot_fitness.sh input.csv output.png

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Default names if not provided
if [ -z "$INPUT_FILE" ]; then
    INPUT_FILE="fitness_log.csv"
fi

if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="fitness_plot.png"
fi

# Check if gnuplot is installed
if ! command -v gnuplot &> /dev/null; then
    echo "Error: gnuplot is not installed."
    exit 1
fi

# Generate gnuplot script
gnuplot -persist <<-EOF
    set datafile separator ","
    set title "Fitness Over Generations"
    set xlabel "Generation"
    set ylabel "Fitness"
    set key outside
    set grid
    set term pngcairo size 1000,600
    set output "${OUTPUT_FILE}"

    plot \
        "${INPUT_FILE}" using 1:2 with linespoints title "Min Fitness" lc rgb "blue", \
        "${INPUT_FILE}" using 1:3 with linespoints title "Max Fitness" lc rgb "red", \
        "${INPUT_FILE}" using 1:4 with linespoints title "Avg Fitness" lc rgb "green"
EOF

echo "Plot saved as ${OUTPUT_FILE}"
