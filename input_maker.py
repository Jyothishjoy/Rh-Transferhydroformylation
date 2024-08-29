import os
import subprocess

# Define the main directory path
main_directory = os.getcwd()  # Replace with your actual path
# Define the path to the template files
in_template_file_path = os.path.join(main_directory, 'inp_template.txt')  # Replace with the actual path

# Function to create subdirectories and process .xyz files
def process_xyz_files(directory, in_template_path):
    with open(in_template_path, 'r') as template:
        template_content = template.read()
    for root, _, files in os.walk(directory):
        for filename in files:
            if filename.endswith('.xyz'):
                xyz_file = os.path.join(root, filename)
                file_name, _ = os.path.splitext(filename)
                # Create a subdirectory with the XYZ file name
                subdirectory = os.path.join(root, file_name)
                os.makedirs(subdirectory, exist_ok=True)
                # Create and write the .inp file in the subdirectory
                in_file = os.path.join(subdirectory, f'{file_name}.com')
                with open(in_file, 'w') as in_template:
                    in_template.write(template_content)
                    with open(xyz_file, 'r') as xyz:
                        lines = xyz.readlines()[2:]
                    in_template.writelines(lines)
                    in_template.write('\n\n\n')
if __name__ == "__main__":
    # process_mol_files(main_directory)  # Uncomment this line if needed
    process_xyz_files(main_directory, in_template_file_path)


