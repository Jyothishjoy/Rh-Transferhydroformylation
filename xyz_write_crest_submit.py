import os
import glob
import subprocess
import shutil

from AaronTools.geometry import Geometry
from AaronTools.fileIO import FileReader
from AaronTools.utils.utils import get_filename, glob_files
from AaronTools.theory import *

current_dir = os.getcwd()
sub_script = os.path.join(current_dir, "submit_crest.sh")

for f in glob.glob("*.log", recursive=True):
    try:
        dirname = f.split(".")[0]
        infile = FileReader(f, just_geom=False)
        geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
        
        output_dir = dirname
        os.makedirs(output_dir, exist_ok=True)
        
        geom.write(outfile=os.path.join(output_dir, f"{dirname}.xyz"))
        shutil.copy(sub_script, output_dir)
        
        # Run the job from the output directory
        subprocess.run("sbatch submit_crest.sh", shell=True, cwd=output_dir)
        
        print(f"Successfully processed {f}")

    except Exception as e:
        print(f"Failed to process {f}: {e}")


