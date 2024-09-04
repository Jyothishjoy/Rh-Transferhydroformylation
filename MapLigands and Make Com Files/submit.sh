#!/bin/bash

# Set the main directory path
main_dir=$(pwd)

# Loop through each subdirectory in the main directory
for subdir in "$main_dir"/*/; do
    # Check if there are any .com files in the subdirectory
    if compgen -G "$subdir"*.com > /dev/null; then
        # Navigate to the subdirectory
        cd "$subdir" || exit
        
        # Execute the command
        run16 12 24 70
        
        # Go back to the main directory
        cd "$main_dir" || exit
    fi
done


