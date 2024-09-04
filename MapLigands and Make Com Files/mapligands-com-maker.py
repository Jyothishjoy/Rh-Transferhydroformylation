import os
import glob
import subprocess

# Function to create subdirectories and process .xyz files
def process_xyz_files(directory, output_file, headder_template_path, footer_template_path):
    with open(headder_template_path, 'r') as headder_template:
        headder_content = headder_template.read()
    with open(footer_template_path, 'r') as footer_template:
        footer_content = footer_template.read()

    for root, _, files in os.walk(directory):
        for filename in files:
            if filename == output_file:
                xyz_file = os.path.join(root, filename)
                file_name, _ = os.path.splitext(filename)
                # Create a subdirectory with the XYZ file name
                subdirectory = os.path.join(root, file_name)
                os.makedirs(subdirectory, exist_ok=True)
                # Create and write the .com file in the subdirectory
                in_file = os.path.join(subdirectory, f'{file_name}.com')
                with open(in_file, 'w') as in_template:
                    in_template.write(headder_content)
                    with open(xyz_file, 'r') as xyz:
                        lines = xyz.readlines()[2:]  # Skip the first two lines
                    in_template.writelines(lines)
                    in_template.write('\n')
                    in_template.write(footer_content)
                    in_template.write('\n\n\n')


# Main directory
main_dir = os.getcwd()
xyz_files = glob.glob(os.path.join(main_dir, "*.xyz"))
ligand = "EZAZEK_Rh_1_2_PP_98"

# Define template paths for different tasks
template_paths = {
    'OxiAdnTS': (os.path.join(main_dir, 'OxiAdnTS_headder.txt'), os.path.join(main_dir, 'OxiAdnTS_footer.txt')),
    'AcylTS': (os.path.join(main_dir, 'AcylTS_headder.txt'), os.path.join(main_dir, 'AcylTS_footer.txt')),
    'betaHTS': (os.path.join(main_dir, 'betaHTS_headder.txt'), os.path.join(main_dir, 'betaHTS_footer.txt')),
    'COextTS': (os.path.join(main_dir, 'COextTS_headder.txt'), os.path.join(main_dir, 'COextTS_footer.txt')),
}

# Logic to select the appropriate template based on some condition (e.g., file name)
for xyz_file in xyz_files:
    output_file = f"{os.path.splitext(os.path.basename(xyz_file))[0]}_{ligand.split('_')[0]}.xyz"
    # Run the mapLigand process
    subprocess.run(["mapLigand.py", xyz_file, "-l", f"1,2={ligand}", "-o", output_file], check=True)

        # Example: Based on file name or other criteria, select the correct template
    if 'OxiAddn-TS' in xyz_file:
        headder_path, footer_path = template_paths['OxiAdnTS']
    elif 'Acyl-TS' in xyz_file:
        headder_path, footer_path = template_paths['AcylTS']
    elif 'b-hydride-Elim-TS' in xyz_file:
        headder_path, footer_path = template_paths['betaHTS']
    elif 'CO-extrusion-TS' in xyz_file:
        headder_path, footer_path = template_paths['COextTS']
    else:
        # Fallback or default case
        headder_path, footer_path = os.path.join(main_dir, 'inp_template_headder.txt'), os.path.join(main_dir, 'inp_template_footer.txt')

    
    # Process the .xyz files
    process_xyz_files(main_dir, output_file, headder_path, footer_path)

