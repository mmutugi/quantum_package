#!/usr/bin/env python3

import os
import subprocess

xyz_files = [f for f in os.listdir('.') if f.endswith('.xyz')]

if not xyz_files:
    print("No .xyz files found in the current directory.")
    exit(1)

for xyz_file in xyz_files:
    base_name = os.path.splitext(xyz_file)[0]
    ezfio_file = f"{base_name}.ezfio"

    try:
        
        subprocess.run(
            f"qp_create_ezfio -b cc-pvdz {xyz_file} -m 1 -o {ezfio_file}",
            shell=True,
            check=True
        )
        
        print(f"SCF completed for {ezfio_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred during EZFIO creation or SCF run for {base_name}: {e}")
        continue


ezfio_files = [f for f in os.listdir('.') if f.endswith('.ezfio')]

for ezfio_file in ezfio_files:
    job_name = ezfio_file[:10] + "_job"

    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=16G
#SBATCH --account=col174
#SBATCH -t 1:00:00

export OMP_NUM_THREADS=128
export OMP_PROC_BIND=false
export QP_MAXMEM=250

source ~/qp2/quantum_package.rc
cd $SLURM_SUBMIT_DIR

qp set perturbation pt2_relative_error 5.e-3
qp set determinants n_det_max 5e6
qp set davidson distributed_davidson False
srun -N 1 -n 1 qp_run scf {ezfio_file} > {ezfio_file}.scf.out
srun -N 1 -n 1 qp_run fci {ezfio_file} > {ezfio_file}.fci.out
"""

    slurm_filename = f"run_{ezfio_file}.slurm"
    with open(slurm_filename, "w") as slurm_file:
        slurm_file.write(slurm_script)

    try:
        print(f"Submitting SLURM job for {ezfio_file}...")
        subprocess.run(
            f"sbatch {slurm_filename}",
            shell=True,
            check=True
        )
        print(f"SLURM job submitted for {ezfio_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while submitting SLURM job for {ezfio_file}: {e}")

print("All SLURM jobs submitted!")

