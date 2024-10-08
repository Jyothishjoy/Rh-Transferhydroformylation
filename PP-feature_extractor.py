import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
import re
from AaronTools.finders import BondedElements, BondedTo, WithinRadiusFromAtom, NotAny
from AaronTools.internal_coordinates import Bond
from AaronTools.geometry import Geometry
from AaronTools.component import Component
from AaronTools.fileIO import FileReader
from AaronTools.finders import NotAny, BondedTo

# Load the data
data = pd.read_excel("Experimental_results.xlsx", sheet_name='Sheet2')
subdirs = data['Ligand']
current_directory = os.getcwd()

# Define the extraction function
def parse_steric(subdir_path):
    ConeAngle = None
    BuriedVol_3 = None
    BuriedVol_4 = None
    BuriedVol_5 = None
    BuriedVol_6 = None
    BuriedVol_7 = None
    B1 = None
    B2 = None
    B3 = None
    B4 = None
    B5 = None
    L = None
    error = None

    log_file_path = os.path.join(subdir_path, f'{subdir}_NBO.log')
    
    if not log_file_path:
        error = "Log file not found"
        return ConeAngle, BuriedVol_3,BuriedVol_4, BuriedVol_5, BuriedVol_6, BuriedVol_7, B1, B2, B3, B4, B5, L, error
    
    try:
        infile = FileReader(log_file_path, just_geom=False)
        geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
        metal = geom.find("Pd")[0]
        CO_carbon = geom.find(BondedElements("Pd", "O"), "C")[0]
        Pd_fragment = geom.get_fragment(metal, stop=CO_carbon, as_object=True)
        
        ligand_atoms = Pd_fragment.find(NotAny(metal))

        # Calculate cone angle
        comp = Component(ligand_atoms)
        key_atoms = comp.find(BondedTo(metal))  # Find coordinating atoms
        for k in key_atoms:
            k.tags.add("key")
        
        cone_angle = comp.cone_angle(center=metal, method='exact', radii='umn')
        ConeAngle = cone_angle  # Assign the calculated cone angle

    except Exception as e:
        #raise e
        error = str(e)
    
    try:
        # Calculate %Buried Volume @3A
        Pd_buried_volume_3 = Pd_fragment.percent_buried_volume(metal, radius=3.0, radii='umn') 
        BuriedVol_3 = Pd_buried_volume_3 
        # Calculate %Buried Volume @ 4A
        Pd_buried_volume_4 = Pd_fragment.percent_buried_volume(metal, radius=4.0, radii='umn') 
        BuriedVol_4 = Pd_buried_volume_4
        # Calculate %Buried Volume @ 5A
        Pd_buried_volume_5 = Pd_fragment.percent_buried_volume(metal, radius=5.0, radii='umn') 
        BuriedVol_5 = Pd_buried_volume_5
        # Calculate %Buried Volume @ 6A
        Pd_buried_volume_6 = Pd_fragment.percent_buried_volume(metal, radius=6.0, radii='umn') 
        BuriedVol_6 = Pd_buried_volume_6
        # Calculate %Buried Volume @ 7A
        Pd_buried_volume_7 = Pd_fragment.percent_buried_volume(metal, radius=7.0, radii='umn') 
        BuriedVol_7 = Pd_buried_volume_7
    
    except Exception as e:
        error = str(e)

    try:
        candidates = geom.find(WithinRadiusFromAtom(metal, 4.0), NotAny(metal))
        coords = geom.coordinates(candidates)

        distances = np.linalg.norm(coords - metal.coords, axis=1)
        d1, d2 = np.partition(distances, 1)[:2]

        # Extract the elemental symbol of the atoms connected to the metal
        atom1 = candidates[np.where(distances == d1)[0][0]]
        atom2 = candidates[np.where(distances == d2)[0][0]]

        CO = geom.find("C", BondedElements("Pd", "O", match_exact=False), [atom1, atom2])[0]
        not_CO = geom.find([atom1, atom2], NotAny(CO))
        L_axis = CO.bond(metal)
        ligand_atoms = set([])
        for atom in not_CO:
            ligand_atoms.update(geom.get_fragment(atom, stop=metal))

        sterimol = geom.sterimol(L_axis=L_axis, start_atom=metal, targets=NotAny(metal), radii="umn")
        B1 = sterimol['B1']
        B2 = sterimol['B2']
        B3 = sterimol['B3']
        B4 = sterimol['B4']
        B5 = sterimol['B5']
        L = sterimol['L']

    except Exception as e:
        error = str(e)

    return ConeAngle, BuriedVol_3, BuriedVol_4, BuriedVol_5, BuriedVol_6, BuriedVol_7, B1, B2, B3, B4, B5, L, error

def parse_structural(subdir_path):
    BiteAngle = None
    Rh_P_min_dist = None
    Rh_P_max_dist = None
    Rh_P_avg_dist = None
    P_P_dist = None
    R_Rh_R_Angle_avg = None
    error = None

    log_file_path = os.path.join(subdir_path, f'{subdir}_NBO.log')
    
    if not log_file_path:
        error = "Log file not found"
        return BiteAngle, Rh_P_min_dist, Rh_P_max_dist, Rh_P_avg_dist,P_P_dist, R_Rh_R_Angle_avg, error
    
    try:
        infile = FileReader(log_file_path, just_geom=False)
        geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
        metal = geom.find("Pd")[0]
        P1 = geom.find("P")[0]
        P2 = geom.find("P")[1]
        
        dist_metal_P1 = metal.dist(P1)
        dist_metal_P2 = metal.dist(P2)
        min_dist = np.min([dist_metal_P1, dist_metal_P2])
        Rh_P_min_dist = min_dist
        max_dist = np.max([dist_metal_P1, dist_metal_P2])
        Rh_P_max_dist = max_dist
        avg_dist = np.mean([dist_metal_P1, dist_metal_P2])
        Rh_P_avg_dist = avg_dist
        PP_dist = P1.dist(P2)
        P_P_dist = PP_dist
        
        angle_radians = geom.angle(P1, metal, P2)
        # Convert radians to degrees and round to 2 digits
        bite_angle_degrees = round(np.degrees(angle_radians), 2)
        BiteAngle = bite_angle_degrees
    
    except Exception as e:
        error = str(e)
        
    try:
        infile = FileReader(log_file_path, just_geom=False)
        geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
        metal = geom.find("Pd")[0]
        CO_carbon = geom.find(BondedElements("Pd", "O"), "C")[0]
        Pd_fragment = geom.get_fragment(metal, stop=CO_carbon, as_object=True)
        P1 = Pd_fragment.find("P")[0]
        P2 = Pd_fragment.find("P")[1]
        
        bonded_to_P1 = Pd_fragment.find(BondedTo(P1), "C")
        angle_radians_P1_1 = geom.angle(bonded_to_P1[0], metal, bonded_to_P1[1])
        angle_degrees_P1_1 = round(np.degrees(angle_radians_P1_1), 2)
        angle_radians_P1_2 = geom.angle(bonded_to_P1[1], metal, bonded_to_P1[2])
        angle_degrees_P1_2 = round(np.degrees(angle_radians_P1_2), 2)
        angle_radians_P1_3 = geom.angle(bonded_to_P1[0], metal, bonded_to_P1[2])
        angle_degrees_P1_3 = round(np.degrees(angle_radians_P1_3), 2)
        avg_P1_angle = round(np.mean([angle_degrees_P1_1, angle_degrees_P1_2, angle_degrees_P1_3]), 2)
          
        bonded_to_P2 = Pd_fragment.find(BondedTo(P2), "C")
        angle_radians_P2_1 = geom.angle(bonded_to_P2[0], metal, bonded_to_P2[1])
        angle_degrees_P2_1 = round(np.degrees(angle_radians_P2_1), 2)
        angle_radians_P2_2 = geom.angle(bonded_to_P2[1], metal, bonded_to_P2[2])
        angle_degrees_P2_2 = round(np.degrees(angle_radians_P2_2), 2)
        angle_radians_P2_3 = geom.angle(bonded_to_P2[0], metal, bonded_to_P2[2])
        angle_degrees_P2_3 = round(np.degrees(angle_radians_P2_3), 2)
        avg_P2_angle = round(np.mean([angle_degrees_P2_1, angle_degrees_P2_2, angle_degrees_P2_3]), 2)
        
    except Exception as e:
        raise e
         
    return BiteAngle, Rh_P_min_dist, Rh_P_max_dist, Rh_P_avg_dist, P_P_dist, avg_P1_angle, avg_P2_angle

def extract_scf_energy_from_log(log_file_path):
    SCF_energy = None
    with open(log_file_path, 'r') as file:
        for line in file:
            if 'SCF Done' in line:
                SCF_energy = float(line.split('=')[-1].split()[0])  # Update last_energy here
    return SCF_energy

def calculate_frozen_bde_from_log(subdir_path):
    complex_path = os.path.join(subdir_path, f'{subdir}_NBO.log')
    frag1_path = os.path.join(subdir_path, f'{subdir}_Lig_NBO.log')
    frag2_path = os.path.join(current_directory, 'Pd_CO', 'Pd_CO_NBO.log')
    
    if os.path.exists(complex_path) and os.path.exists(frag1_path) and os.path.exists(frag2_path):
        complex_energy = extract_scf_energy_from_log(complex_path)
        frag1_energy = extract_scf_energy_from_log(frag1_path)
        frag2_energy = extract_scf_energy_from_log(frag2_path)

        if complex_energy is not None and frag1_energy is not None and frag2_energy is not None:
            BDE = (frag1_energy + frag2_energy - complex_energy) * 627.503
        return BDE
    else:
        return None

def extract_CO_freq(subdir_path):
    log_file_path = os.path.join(subdir_path, f'{subdir}.log')
    infile = FileReader(log_file_path, just_geom=False)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
    metal = geom.find("Pd")[0]
    CO = geom.find(BondedElements("Pd", "O"), "C")[0]
    bond = Bond(geom.atoms.index(metal), geom.atoms.index(CO))
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

    return vib_freq

def extract_somo_lumo(subdir_path):
    log_file_path = os.path.join(subdir_path, f'{subdir}_NBO.log')
    
    if not log_file_path:
        return None, None, None  

    with open(log_file_path, 'r') as f:
        lines = f.readlines()

        last_alpha_line = None
        homo = None
        lumo = None
        homo_lumo = None
        
        try:
            for line in lines:
                if "Alpha  occ. eigenvalues" in line:
                    last_alpha_line = line
                    continue

                if last_alpha_line and "Alpha virt. eigenvalues" in line:
                    lumo = float(line.split()[4])
                    homo = float(last_alpha_line.split()[-1])
                    homo_lumo = lumo - homo
                    break
        except Exception as e:
            error = str(e)
            print(f"An error occurred: {error}")
            return None, None, None  # Return None values if an error occurs

    return homo, lumo, homo_lumo

def extract_dipole(subdir_path):
    log_file_path = os.path.join(subdir_path, f'{subdir}_NBO.log')
    
    with open(log_file_path, 'r') as file:
        content = file.read()

        start_marker = "Dipole moment (field-independent basis, Debye):"
        end_marker = "Quadrupole moment (field-independent basis, Debye-Ang):"

        start_index = content.rfind(start_marker)
        end_index = content.find(end_marker, start_index)

        if start_index != -1 and end_index != -1:
            section = content[start_index:end_index]

            dipole_line = None
            for line in section.split('\n'):
                if 'Tot=' in line:
                    dipole_line = line
                    break
            if dipole_line:
                items = dipole_line.split()
                dipole = float(items[-1])
            else:
                dipole = None
        else:
            dipole = None

        return dipole

def extract_P_NMR(subdir_path):
    log_file_path = os.path.join(subdir_path, f'{subdir}_NBO.log')
    
    infile = FileReader(log_file_path)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
    metal = geom.find("Pd")[0]
    P1 = geom.find("P")[0]
    P2 = geom.find("P")[1]
    P1_id = geom.atoms.index(P1) +1
    P2_id = geom.atoms.index(P2) +1
    
    with open(log_file_path, 'r') as file:
        content = file.read()

        start_marker = "SCF GIAO Magnetic shielding tensor (ppm):"
        end_marker = " **********************************************************************"

        start_index = content.rfind(start_marker)
        end_index = content.find(end_marker, start_index)

        if start_index != -1 and end_index != -1:
            section = content[start_index:end_index]

            P1_nmr_line = None
            P2_nmr_line = None
            for line in section.split('\n'):
                if f'{P1_id}  P' in line:
                    P1_nmr_line = line
                if f'{P2_id}  P' in line:
                    P2_nmr_line = line
                    break
            if P1_nmr_line:
                items = P1_nmr_line.split()
                P1_Isotropic = float(items[4])
                P1_Anisotropic = float(items[-1])
            if P2_nmr_line:
                items = P2_nmr_line.split()
                P2_Isotropic = float(items[4])
                P2_Anisotropic = float(items[-1])
            else:
                P1_Isotropic, P1_Anisotropic, P2_Isotropic, P2_Anisotropic = None, None, None, None
        else:
            P1_Isotropic, P1_Anisotropic, P2_Isotropic, P2_Anisotropic = None, None, None, None

        return P1_Isotropic, P1_Anisotropic, P2_Isotropic, P2_Anisotropic


def extract_NBO_pop(subdir_path):
    log_file_path = os.path.join(subdir_path, f'{subdir}_Lig_NBO.log')
    
    # Check if log file exists
    if not log_file_path:
        print(f"No log file found for {subdir}.")
        return None, None, None, None

    infile = FileReader(log_file_path)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
    P1 = geom.find("P")[0]
    P2 = geom.find("P")[1]
    P1_id = geom.atoms.index(P1) + 1
    P2_id = geom.atoms.index(P2) + 1
    
    with open(log_file_path, 'r') as file:
        content = file.read()

        start_marker = "Summary of Natural Population Analysis:"
        end_marker = " ======================================================================="

        start_index = content.rfind(start_marker)
        end_index = content.find(end_marker, start_index)

        if start_index != -1 and end_index != -1:
            section = content[start_index:end_index]

            P1_line = None
            P2_line = None
            
            # Create regex patterns for P1 and P2 IDs
            P1_pattern = re.compile(rf'\bP\s+{P1_id}\b')  # \b ensures it matches whole number
            P2_pattern = re.compile(rf'\bP\s+{P2_id}\b')
            
            for line in section.split('\n'):
                # Check for P1
                if P1_pattern.search(line):
                    P1_line = line.strip()
                    break
                
            for line in section.split('\n'):
                # Check for P1
                if P2_pattern.search(line):
                    P2_line = line.strip()
                    break

            if P1_line:
                P1_items = P1_line.split()
                P1_natural_charge = float(P1_items[2])
                P1_natural_population = float(P1_items[-1])
            else:
                P1_natural_charge = None
                P1_natural_population = None

            if P2_line:
                P2_items = P2_line.split()
                P2_natural_charge = float(P2_items[2])
                P2_natural_population = float(P2_items[-1])
            else:
                P2_natural_charge = None
                P2_natural_population = None
        else:
            P1_natural_charge = None
            P1_natural_population = None
            P2_natural_charge = None
            P2_natural_population = None

    return P1_natural_charge, P1_natural_population, P2_natural_charge, P2_natural_population


def extract_NBO_occ(subdir_path):
    log_file_path = os.path.join(subdir_path, f'{subdir}_Lig_NBO.log')
    
    # Check if log file exists
    if not log_file_path:
        print(f"No log file found for {subdir}.")
        return None, None

    infile = FileReader(log_file_path)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
    P1 = geom.find("P")[0]
    P2 = geom.find("P")[1]
    P1_id = geom.atoms.index(P1) + 1
    P2_id = geom.atoms.index(P2) + 1
    
    with open(log_file_path, 'r') as file:
        content = file.read()

        start_marker = "Natural Bond Orbitals (Summary):"
        end_marker = "       -------------------------------"

        start_index = content.rfind(start_marker)
        end_index = content.find(end_marker, start_index)

        if start_index != -1 and end_index != -1:
            section = content[start_index:end_index]

            P1_line = None
            P2_line = None

            # Create regex patterns for P1 and P2 IDs
            P1_pattern = re.compile(rf'LP.*P\s+{P1_id}\b')
            P2_pattern = re.compile(rf'LP.*P\s+{P2_id}\b')
            
            # Search through the section for both P1 and P2 lines
            for line in section.split('\n'):
                if P1_pattern.search(line):
                    P1_line = line.strip()
                if P2_pattern.search(line):
                    P2_line = line.strip()
                if P1_line and P2_line:
                    break

            # Extract lone pair occupancy data
            if P1_line:
                P1_items = P1_line.split()
                P1_LP_Occ = float(P1_items[6])
            else:
                P1_LP_Occ = None

            if P2_line:
                P2_items = P2_line.split()
                P2_LP_Occ = float(P2_items[6])
            else:
                P2_LP_Occ = None

            # Compute average occupancy
            avg_P_LP_occ = np.mean([P1_LP_Occ, P2_LP_Occ])

        else:
            P1_LP_Occ = None
            P2_LP_Occ = None
            avg_P_LP_occ = None
        
    return avg_P_LP_occ
   
    
def parse_electronic(subdir_path):
    BDE = calculate_frozen_bde_from_log(subdir_path)
    TEP_uncorrected = extract_CO_freq(subdir_path)
    HOMO, LUMO, HOMO_LUMO = extract_somo_lumo(subdir_path)
    dipole = extract_dipole(subdir_path)
    P1_Isotropic, P1_Anisotropic, P2_Isotropic, P2_Anisotropic = extract_P_NMR(subdir_path)
    P1_natural_charge, P1_natural_population, P2_natural_charge, P2_natural_population = extract_NBO_pop(subdir_path)
    avg_P_LP_occ = extract_NBO_occ(subdir_path)
    
    return BDE, TEP_uncorrected, HOMO, LUMO, HOMO_LUMO, dipole, P1_Isotropic, P1_Anisotropic, P2_Isotropic, P2_Anisotropic, P1_natural_charge, P1_natural_population, P2_natural_charge, P2_natural_population, avg_P_LP_occ


def extract_conceptual_dft(subdir_path): # https://chemtools.org/sci_doc_conceptual.html
    frag1_path = os.path.join(subdir_path, f'{subdir}_Lig_NBO.log')
    frag1_cat_path = os.path.join(subdir_path, f'{subdir}_Lig_Cat_NBO.log')
    frag1_ani_path = os.path.join(subdir_path, f'{subdir}_Lig_Ani_NBO.log')
    
    if os.path.exists(frag1_path) and os.path.exists(frag1_cat_path) and os.path.exists(frag1_ani_path):
        frag1_energy = extract_scf_energy_from_log(frag1_path)
        frag1_cat_energy = extract_scf_energy_from_log(frag1_cat_path)
        frag1_ani_energy = extract_scf_energy_from_log(frag1_ani_path)

        if frag1_energy is not None and frag1_cat_energy is not None and frag1_ani_energy is not None:
            I = (frag1_cat_energy - frag1_energy)  # Ionization Energy
            A = (frag1_energy - frag1_ani_energy)   # Electron Affinity 
            mue = -(I + A) / 2  # Chemical Potential
            eta = I - A         # Chemical Hardness
            S = 1 / (2 * eta)   # Chemical Softness
            omega = mue**2 / (2 * eta)  # Electrophilicity Index
            omega_negative = (3 * I + A)**2 / (16 * (I - A))
            omega_positive = (I + 3 * A)**2 / (16 * (I - A))

            return I, A, mue, eta, S, omega, omega_negative, omega_positive
        else:
            return None, None, None, None, None, None, None, None


def extract_Fukui_Fun(subdir_path):
    def extract_NBO_POP(file):
        infile = FileReader(file)
        geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
        P1 = geom.find("P")[0]
        P2 = geom.find("P")[1]
        P1_id = geom.atoms.index(P1) + 1
        P2_id = geom.atoms.index(P2) + 1
        
        with open(file, 'r') as f:
            file_content = f.read()

            start_marker = "Summary of Natural Population Analysis:"
            end_marker = "======================================================================="

            start_index = file_content.rfind(start_marker)
            end_index = file_content.find(end_marker, start_index)
        
            if start_index != -1 and end_index != -1:
                section = file_content[start_index:end_index].splitlines()

                P1_line = None
                P2_line = None
            
                # Create regex patterns for P1 and P2 IDs
                P1_pattern = re.compile(rf'\bP\s+{P1_id}\b')
                P2_pattern = re.compile(rf'\bP\s+{P2_id}\b')

                for line in section:
                    if P1_pattern.search(line):
                        P1_line = line.strip()
                        break
                
                for line in section:
                    if P2_pattern.search(line):
                        P2_line = line.strip()
                        break

                P1_natural_population = float(P1_line.split()[-1]) if P1_line else None
                P2_natural_population = float(P2_line.split()[-1]) if P2_line else None

                if P1_natural_population is not None and P2_natural_population is not None:
                    avg_P_pop = np.mean([P1_natural_population, P2_natural_population])
                else:
                    avg_P_pop = None
            else:
                print(f"Failed to extract NBO population data from {file}")
                avg_P_pop = None

            return avg_P_pop

    frag1_path = os.path.join(subdir_path, f'{subdir}_Lig_NBO.log')
    frag1_cat_path = os.path.join(subdir_path, f'{subdir}_Lig_Cat_NBO.log')
    frag1_ani_path = os.path.join(subdir_path, f'{subdir}_Lig_Ani_NBO.log')
    
    if os.path.isfile(frag1_path) and os.path.isfile(frag1_cat_path) and os.path.isfile(frag1_ani_path):
        frag1_P_avg_pop = extract_NBO_POP(frag1_path)
        frag1_cat_P_avg_pop = extract_NBO_POP(frag1_cat_path)
        frag1_ani_P_avg_pop = extract_NBO_POP(frag1_ani_path)
        
        if frag1_P_avg_pop is not None and frag1_cat_P_avg_pop is not None and frag1_ani_P_avg_pop is not None:
            Electrophilicity = frag1_ani_P_avg_pop - frag1_P_avg_pop
            Nucleophilicity = frag1_P_avg_pop - frag1_cat_P_avg_pop
            Radical_attack_susceptibility = (frag1_ani_P_avg_pop - frag1_cat_P_avg_pop) * 0.5
        else:
            Electrophilicity, Nucleophilicity, Radical_attack_susceptibility = None, None, None

        return Electrophilicity, Nucleophilicity, Radical_attack_susceptibility
    else:
        print(f"One or more log files are missing for {subdir}")
        return None, None, None


    
               
# Process each ligand and append the results to the dataframe
for subdir in tqdm(subdirs):
    subdir_path = os.path.join(current_directory, subdir)
    
    if os.path.isdir(subdir_path):
    
        print(f"Start Calculation for {subdir}")
        ConeAngle, BuriedVol_3, BuriedVol_4, BuriedVol_5, BuriedVol_6, BuriedVol_7, B1, B2, B3, B4, B5, L, error = parse_steric(subdir_path)
        BiteAngle, Rh_P_min_dist, Rh_P_max_dist, Rh_P_avg_dist, P_P_dist, avg_P1_angle, avg_P2_angle  = parse_structural(subdir_path)
        BDE, TEP_uncorrected, HOMO, LUMO, HOMO_LUMO, dipole, P1_Isotropic, P1_Anisotropic, P2_Isotropic, P2_Anisotropic, P1_natural_charge, P1_natural_population, P2_natural_charge, P2_natural_population, avg_P_LP_occ = parse_electronic(subdir_path)
        I, A, mue, eta, S, omega, omega_negative, omega_positive = extract_conceptual_dft(subdir_path)
        Electrophilicity, Nucleophilicity, Radical_attack_susceptibility = extract_Fukui_Fun(subdir_path)
    

        # Update the corresponding row in the DataFrame
        data.loc[data['Ligand'] == subdir, 'Cone_Angle'] = ConeAngle
        data.loc[data['Ligand'] == subdir, '%Buried Volume_3A'] = BuriedVol_3
        data.loc[data['Ligand'] == subdir, '%Buried Volume_4A'] = BuriedVol_4
        data.loc[data['Ligand'] == subdir, '%Buried Volume_5A'] = BuriedVol_5
        data.loc[data['Ligand'] == subdir, '%Buried Volume_6A'] = BuriedVol_6
        data.loc[data['Ligand'] == subdir, '%Buried Volume_7A'] = BuriedVol_7
        data.loc[data['Ligand'] == subdir, 'B1'] = B1
        data.loc[data['Ligand'] == subdir, 'B2'] = B2
        data.loc[data['Ligand'] == subdir, 'B3'] = B3
        data.loc[data['Ligand'] == subdir, 'B4'] = B4
        data.loc[data['Ligand'] == subdir, 'B5'] = B5
        data.loc[data['Ligand'] == subdir, 'L'] = L
        data.loc[data['Ligand'] == subdir, 'Bite Angle'] = BiteAngle
        data.loc[data['Ligand'] == subdir, 'Rh-P_min_dist'] = Rh_P_min_dist
        data.loc[data['Ligand'] == subdir, 'Rh-P_max_dist'] = Rh_P_max_dist
        data.loc[data['Ligand'] == subdir, 'Rh-P_avg_dist'] = Rh_P_avg_dist
        data.loc[data['Ligand'] == subdir, 'P-P_dist'] = P_P_dist
        data.loc[data['Ligand'] == subdir, 'Avg C_angle around P1'] = avg_P1_angle
        data.loc[data['Ligand'] == subdir, 'Avg C_angle around P2'] = avg_P2_angle
        data.loc[data['Ligand'] == subdir, 'Metal-Lig-BDE'] = BDE
        data.loc[data['Ligand'] == subdir, 'TEP(uncorrected)'] = TEP_uncorrected
        data.loc[data['Ligand'] == subdir, 'HOMO'] = HOMO
        data.loc[data['Ligand'] == subdir, 'LUMO'] = LUMO
        data.loc[data['Ligand'] == subdir, 'HOMO-LUMO'] = HOMO_LUMO
        data.loc[data['Ligand'] == subdir, 'Complex_Dipole'] = dipole
        data.loc[data['Ligand'] == subdir, 'P1_NMR_Isotropic'] = P1_Isotropic
        data.loc[data['Ligand'] == subdir, 'P1_NMR_Anisotropic'] = P1_Anisotropic
        data.loc[data['Ligand'] == subdir, 'P2_NMR_Isotropic'] = P2_Isotropic
        data.loc[data['Ligand'] == subdir, 'P2_NMR_Anisotropic'] = P2_Anisotropic
        data.loc[data['Ligand'] == subdir, 'P1_natural_charge'] = P1_natural_charge
        data.loc[data['Ligand'] == subdir, 'P1_natural_population'] = P1_natural_population
        data.loc[data['Ligand'] == subdir, 'P2_natural_charge'] = P2_natural_charge
        data.loc[data['Ligand'] == subdir, 'P2_natural_population'] = P2_natural_population
        data.loc[data['Ligand'] == subdir, 'P_LP_Occ_avg'] = avg_P_LP_occ
        data.loc[data['Ligand'] == subdir, 'Ionization Energy'] = I
        data.loc[data['Ligand'] == subdir, 'Electron Affinity'] = A
        data.loc[data['Ligand'] == subdir, 'Chemical Potential'] = mue
        data.loc[data['Ligand'] == subdir, 'Chemical Hardness'] = eta
        data.loc[data['Ligand'] == subdir, 'Chemical Softness'] = S
        data.loc[data['Ligand'] == subdir, 'Electrophilicity'] = Electrophilicity
        data.loc[data['Ligand'] == subdir, 'Nucleophilicity'] = Nucleophilicity
        data.loc[data['Ligand'] == subdir, 'Radical_attack_susceptibility'] = Radical_attack_susceptibility


        # Save the updated dataframe back to the CSV file after each ligand is processed
        data.to_csv("extracted_data.csv", index=False)

        # Print errors if they occur
        if error:
            print(f"Error with {subdir}: {error}")

    print(f"Processing completed. Results saved to extracted_data.csv.")
