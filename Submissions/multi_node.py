#!/usr/bin/env python3
import os
import subprocess
import numpy as np

ezfio_files = [f for f in os.listdir('.') if f.endswith('.ezfio')]

for ezfio_file in ezfio_files:
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={ezfio_file[:10]}_job  
#SBATCH --partition=compute
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --account=col174
#SBATCH -t 24:00:00

export OMP_NUM_THREADS=128
export OMP_PROC_BIND=false
export QP_MAXMEM=250

source ~/qp2/quantum_package.rc

cd $SLURM_SUBMIT_DIR
qp set_file {ezfio_file}
qp_run save_natorb {ezfio_file} > {ezfio_file}.natorbs.out
qp set determinants s2_eig False
qp set determinants read_wf True
qp set determinants n_det_max 5e8
qp set davidson distributed_davidson True
qp_srun fci {ezfio_file} > {ezfio_file}.fci.natorbs.out
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


