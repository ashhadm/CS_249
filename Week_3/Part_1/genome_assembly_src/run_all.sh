#!/bin/bash
# run_all.sh - Run all parts of Task 1.3

# Run setup
echo "Running setup..."
bash setup.sh

# Run each part
echo -e "\n\n========== Running Task 1.3 Part 1 ==========\n"
python3 task1_3_part1.py

echo -e "\n\n========== Running Task 1.3 Part 2 ==========\n"
python3 task1_3_part2.py

echo -e "\n\n========== Running Task 1.3 Part 3 ==========\n"
python3 task1_3_part3.py

echo -e "\n\n========== Running Task 1.3 Part 4 ==========\n"
python3 task1_3_part4.py

echo -e "\n\nAll tasks completed!"
echo "Results are stored in the genome_assembly_output directory."
