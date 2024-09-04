# How to use AaronTool's MapLigand Functionality to do MultiSubstitutions

mapligands.py script uses AaronTool's MapLigand functionality to mount all bidentate ligands from AaronTools library into the PP-template.xyz file. 

## How to add new ligands to the Aaron_library?
https://aarontools.readthedocs.io/en/latest/other_docs/libraries.html

Windows:
One can add desired ligands to Aaron tools custom library located under "C:\Users\jyothish\Aaron_libs"

There need to have some preprocessing before adding ligands to the library.

1. If you are starting from mol files (eg: in ReaLigands), convert them first into xyz files using openbabel with the following command;

obabel *.mol -oxyz -m


2. Thus generated XYZ file retains (luckly!) the ligand connectivity info from ReaLigands as the title of the XYZ file. 
This title info can be reformatted as an Aaron tools recognizable data. Just reformat the connectivity info in the following format. K:1,2; where, 1 and 2 are the connecting atoms. 
Use python to automate this.

import os
directory_path = os.getcwd()
for filename in os.listdir(directory_path):
    if filename.endswith('.xyz'):
        file_path = os.path.join(directory_path, filename)
        # Read the content of the file
        with open(file_path, 'r') as file:
            lines = file.readlines()
        # Check if the file has at least two lines
        if len(lines) >= 2:
            # Modify the second line
            lines[1] = f'K:{lines[1].strip()};\n'
            # Write the modified content back to the file
            with open(file_path, 'w') as file:
                file.writelines(lines)
print("Content in the second line of each .xyz file has been replaced.")


3. Now, the correctly formatted xyz files of the ligands can be deposited to the Aaron tools custom library.

4. One can use command line scripts to generate bulk ligand substitutions. Use the following shell script:


#!/bin/env bash

for lig in `mapLigand.py -ls elements:P,P`; do
    mapLigand.py PP-template.xyz -l 1,2=$lig -o $lig.xyz;
done

Save this as a sh file and execute. -ls lists all the ligands, and sort based on the connecting elements P,P. 
We can also use denticity:2 to sort all the bidentate ligands. 
Then execute mapLigands on top of the selected ones using mapLigand.py template.xyz -l 1,2=$lig -o $lig.xyz. 
-l represents the ligand, 1,2 are the positions for the bidentate ligand to get replaced at the template file.

5. This script will generate many complexes in the xyz format.
input_maker.py will use the inp_template.txt and above generated xyz files to create Gaussian input files.

6. submit.sh script will submit the Gaussian .com files to the supercomputer.
