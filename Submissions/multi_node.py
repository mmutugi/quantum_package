#!/usr/bin/env python3
import os
import subprocess
import numpy as np

ezfio_files = [f for f in os.listdir('.') if f.endswith('.ezfio')]

for ezfio_file in ezfio_files:
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={ezfio_file[:10]}_job  # Ensure job name isn't too long
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
source ~/qp2/quantum_package.rc
cd $SLURM_SUBMIT_DIR
qp set_file {ezfio_file}
qp_run save_natorb {ezfio_file}
qp set determinants read_wf true
qp set determinants n_det_max 10000000
qp_srun fci {ezfio_file} > {ezfio_file}_natorbs.out
"""

    slurm_filename = f"run_{ezfio_file}.slurm"
    with open(slurm_filename, "w") as slurm_file:
        slurm_file.write(slurm_script)

    try:
        print(f"Submitting SLURM job for {ezfio_file}...")
        subprocess.run(
            f"sbatch {slurm_filename}",
            shell=True, check=True
        )
        print(f"SLURM job submitted for {ezfio_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while submitting SLURM job for {ezfio_file}: {e}")

print("All SLURM jobs submitted!")

