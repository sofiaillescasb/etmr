#!/usr/bin/env bash

# Exit on error
set -euo pipefail

INPUT_DIR="/home/sofia/Projects/etmr/defaria/snRNA/data/processed/merged"
PYTHON_SCRIPT="03_integrate.py"

# Check that python is available (from active conda env)
if ! command -v python &> /dev/null; then
    echo "Error: python not found in PATH. Activate your conda environment first."
    exit 1
fi

# Optional: print which python is being used
echo "Using python: $(which python)"

# Loop through each file in the directory
for file in "${INPUT_DIR}"/*; do
    # Skip if not a regular file
    [[ -f "$file" ]] || continue

    echo "Processing: $file"
    python "$PYTHON_SCRIPT" "$file"
done

echo "All files processed."
