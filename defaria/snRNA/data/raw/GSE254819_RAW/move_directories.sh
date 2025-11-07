#!/bin/bash

# Loop through all files ending with _barcodes.tsv.gz
for file in *barcodes.tsv.gz; do
  # Extract prefix before _barcodes.tsv.gz
  prefix="${file%%_barcodes.tsv.gz}"

  # Create the directory (if not already exists)
  mkdir -p "$prefix"

# Move all files starting with that prefix into the directory
  for f in "${prefix}"*; do
    # Remove the prefix from the filename
    newname="${f#$prefix}"
    # Remove leading underscore if present
    newname="${newname#_}"
    # Move and rename the file
    mv "$f" "$prefix/$newname"
    echo "Moved and renamed '$f' -> '$prefix/$newname'"
    done

done
