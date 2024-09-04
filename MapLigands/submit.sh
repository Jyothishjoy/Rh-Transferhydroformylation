#!/bin/bash

# Set the main directory path
main_dir="/home/fslcollab286/nobackup/archive/Retrohydroformylation/Xantphos-Co-Rh-Ir-Daniel/Tolman_Electronic_parameters/B3LYP_D3BJ_Def2SVP/ReaLigands_TEP"

# Loop through each subdirectory in the main directory
for subdir in "$main_dir"/*/; do
    # Check if there are any .com files in the subdirectory
    if compgen -G "$subdir"*.com > /dev/null; then
        # Navigate to the subdirectory
        cd "$subdir" || exit
        
        # Execute the command
        sbatch *.sh
        
        # Go back to the main directory
        cd "$main_dir" || exit
    fi
done


