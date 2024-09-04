# Functionality of mapligands-com-maker.py Script

`mapligands-com-maker.py` script can automate the preparation of all the input files required for the calculation of the transfer hydroformylation mechanism for a given PP-bidentate ligand. 

This script uses the core structures of every point in the energy profile. Then use AaronTool's `mapLigand` functionality to mound PP-ligand of our choice to the core template.
In the second step, the generated XYZ files are taken, a header and footer are added to it, and saved as a .com file for G16 calculation.
