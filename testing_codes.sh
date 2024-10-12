#!/bin/bash

# Test whether all codes run via julia terminal

# Directory containing the Julia files
DIRECTORY="code"

# Check if the directory exists
if [ ! -d "$DIRECTORY" ]; then
  echo "Directory $DIRECTORY does not exist."
  exit 1
fi

# Find all .jl files in the directory and execute them with Julia
for file in "$DIRECTORY"/*.jl; do
  if [ -f "$file" ]; then
    echo "Running $file"
    julia "$file"
  else
    echo "No Julia files found in $DIRECTORY."
    exit 1
  fi
done