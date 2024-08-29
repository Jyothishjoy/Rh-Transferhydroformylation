import os
import csv
import numpy as np
from AaronTools.geometry import Geometry
from AaronTools.fileIO import FileReader
from AaronTools.finders import BondedElements
from AaronTools.internal_coordinates import Bond
from AaronTools.utils.utils import get_filename
import glob

def extract_CO_freq(log_file_path):
    infile = FileReader(log_file_path, just_geom=False)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
    metal = geom.find("Pd")[0]
    CO = geom.find(BondedElements("Pd", "O"), "C")[0]
    bond = Bond(geom.atoms.index(metal), geom.atoms.index(CO))
    distance = CO.dist(metal)
    s_vec = bond.s_vector(geom.coords)
    freq = infile["frequency"]

    max_overlap = None
    best_stretch = None
    for mode in freq.data:
        overlap = abs(np.dot(s_vec, np.reshape(mode.vector, -1)))
        if max_overlap is None or overlap > max_overlap:
            max_overlap = overlap
            best_stretch = mode

    vib_freq = best_stretch.frequency
    force_const = best_stretch.forcek

    return distance, vib_freq, force_const

def extract_bite_angle(log_file_path):
    infile = FileReader(log_file_path, just_geom=False)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
    metal = geom.find("Pd")[0]
    P1 = geom.find("P")[0]
    P2 = geom.find("P")[1]
    angle_radians = geom.angle(P1, metal, P2)
    # Convert radians to degrees and round to 2 digits
    bite_angle_degrees = round(np.degrees(angle_radians), 2)

    return bite_angle_degrees

# Main function
if __name__ == "__main__":
    # Process each subdirectory and write to CSV
    with open('extracted_TEP_ConeAngl.csv', 'w', newline='') as csvfile:
        fieldnames = [
            'CSD_ID', 'CO_Stretch(UnCorr)', 'Pd-C Bond Length', 'Pd-C Force Constant', 'P_P Bite Angle'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        main_directory = os.getcwd()
        for subdir in os.listdir(main_directory):
            subdir_path = os.path.join(main_directory, subdir)
            if os.path.isdir(subdir_path):
                log_files = glob.glob(os.path.join(subdir_path, '*.log'))
                if log_files:
                    log_file_path = log_files[0]  # Assuming there is only one log file per subdir
                    try:
                        # Extract various data
                        distance, vib_freq, force_const = extract_CO_freq(log_file_path)
                        bite_angle_degrees = extract_bite_angle(log_file_path)
                        
                        # Write the row to the CSV file
                        writer.writerow({
                            'CSD_ID': subdir,
                            'CO_Stretch(UnCorr)': vib_freq,
                            'Pd-C Bond Length' : distance,
                            'Pd-C Force Constant': force_const,
                            'P_P Bite Angle': bite_angle_degrees,
                        })
                        print(f"Processed {subdir}")
                    except Exception as e:
                        print(f"Error processing {subdir}: {e}")

