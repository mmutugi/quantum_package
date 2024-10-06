#!/usr/bin/env python3

import numpy as np
import sys
from scipy.io import FortranFile
from ezfio import ezfio

def read_ezfio(ezfio_filename):
    ezfio.set_file(ezfio_filename)
    
    output_file = 'CI_coeff.dat'
    MCSCF_dat  = 'MCSCF_MO.dat'
    afqmc_in_filename = 'afqmc.in'
    filename = 'V2b_AO_cholesky.mat'
    psi_coeff = ezfio.get_determinants_psi_coef()[0]  
    psi_det = ezfio.get_determinants_psi_det()
    mo_coef = ezfio.get_mo_basis_mo_coef()
    alpha_num = ezfio.get_electrons_elec_alpha_num()
    beta_num = ezfio.get_electrons_elec_beta_num()
    n_basis = ezfio.get_mo_basis_mo_num()
    NORB = n_basis 
    n_electron = beta_num + alpha_num
    n_det  = ezfio.get_determinants_n_det()

    alpha_beta_det = []

    # Loop
    for i, det_pair in enumerate(psi_det):
        try:
            alpha_det = det_pair[0]  
            beta_det = det_pair[1]   
            psi_coef = psi_coeff[i]  

            alpha_det_binary = ''.join([bin(det) for det in alpha_det])  
            beta_det_binary = ''.join([bin(det) for det in beta_det])    

            # Append data
            alpha_beta_det.append([alpha_det_binary, beta_det_binary, psi_coef])

        except IndexError:
            print(f"Warning: Index error at i={i}, not aligned.")
            break

    with open(output_file, 'w') as f:
        f.write(f"{'Alpha det':20} {'Beta det':20} {'Psi_coefficients':20}\n")
        
        for alpha_det, beta_det, psi_coef in alpha_beta_det:
            f.write(f"{alpha_det:20} {beta_det:20} {psi_coef:20.16f}\n")

    print(f"Success writing in {output_file}")

    with open(MCSCF_dat, 'w') as f:
        mo_s = [coef for sublist in mo_coef for coef in sublist]
        f.write(f"MO coefficients: {len(mo_s)}\n")
        for coef in mo_s:
            f.write(f"{coef:24.16f}\n")
        f.write('\n')
                                    
    print("MCSF_MO data also written")

    return NORB, n_basis, n_electron, alpha_num, beta_num, n_det, afqmc_in_filename,filename


def dump_fci(fcidump_file, NORB): 
    eri = np.zeros([NORB, NORB, NORB, NORB])

    with open(fcidump_file, 'r') as file:
        for line in file:
            if not line.strip() or line.startswith('&FCI') or line.startswith(' ORBSYM'):
                continue
            
            parts = line.split()
            
            if len(parts) == 5:
                val = float(parts[0]) 
                i = int(parts[1]) - 1  
                j = int(parts[2]) - 1
                k = int(parts[3]) - 1
                l = int(parts[4]) - 1

                eri[i, j, k, l] = val
    return eri

def core_mcd(eri, NORB, chmax=2000, tolcd=1e-5):
    L = np.zeros((NORB * NORB, chmax))
    Dmax = np.zeros((NORB * NORB))
    ng = 0

    V = eri.reshape(NORB**2, NORB**2) / 2
    Dmax[:] = np.diagonal(V)

    while True:
        nu_max = np.argmax(Dmax)
        vmax = Dmax[nu_max]

        if vmax < tolcd:
            break
        if ng >= chmax:
            print("WARNING -- EXCEEDED CHOLESKY LIMIT -- CHECK CONVERGENCE")
            break

        L[:, ng][:] = V[:, nu_max]

        if ng > 0:
            L[:, ng] -= np.dot(L[:, 0:ng], (L.T)[0:ng, nu_max])
        L[:, ng] /= np.sqrt(vmax)

        for mu in range(NORB * NORB):
            Dmax[mu] -= L[mu, ng] ** 2

        ng += 1
    print("Cholesky fields:", ng)
    return L[:, :ng]

def write_mcd_core(filename, NORB, mcd):
    with FortranFile(filename, 'w') as g:
        g.write_record(2)
        g.write_record(np.array([NORB * (NORB + 1) // 2, mcd.shape[1]], dtype=np.int32))
        Lmat = np.zeros((NORB, NORB))
        Llst = np.zeros((NORB * (NORB + 1) // 2))

        for gamma in range(mcd.shape[1]):
            Lmat = mcd[:, gamma].reshape((NORB, NORB))
            counter = 0
            for ib in range(NORB):
                for jb in range(ib, NORB):
                    Llst[counter] = Lmat[ib, jb]
                    counter += 1
            g.write_record(Llst[:])
        return filename

def afqmc_in(n_basis, n_electron, alpha_num, beta_num, n_det, afqmc_in_filename,filename):
    with open(afqmc_in_filename, 'w') as f:
        #:::::::::::::::::::::::::::::::::

        #:::::::::::::::::::::::::::::::::
        f.write(f"CHEM_SYS \"{'defaults_single_det'}\"\n")
        f.write(f"FLAG_BEG_FP {'false'}\n")
        f.write(f"FLAG_CR_FP {'false'}\n")
        f.write(f"BLK_START_FP  {'5'}\n")
        f.write(f"FLAG_CHOLESKY {'true'}\n")
        f.write(f"FLAG_INIT_RHF_WLK {'false'}\n")
        f.write(f"FLAG_INIT_CAS {'true'}\n")
        f.write(f"FLAG_CHANGE_REF {'false'}\n")
        f.write(f"RAND_SEED  {'0'}\n")
        f.write(f"M_BASIS {n_basis}\n")
        f.write(f"N_ELEC {n_electron}\n")
        f.write(f"ENERGY_FC {'0.0'}\n")
        f.write(f"N_UP {alpha_num}\n")
        f.write(f"N_DN {beta_num}\n")
        f.write(f"N_ACT_ORBS {n_basis}\n")
        f.write(f"N_ACT_UP {alpha_num}\n")
        f.write(f"N_ACT_DN {beta_num}\n")
        f.write(f"N_DET {n_det}\n")

        f.write(f"N_WLK {'4000'}\n")
        f.write(f"N_BLK {'1000'}\n")
        f.write(f"N_BLKSTEPS {'20'}\n")
        f.write(f"ITV_MODSVD {'2'}\n")
        f.write(f"ITV_PC {'20'}\n")
        f.write(f"ITV_EM {'20'}\n")
        
        f.write(f"DELTAU {'0.005'}\n")
        f.write(f"DELTA_EH {'20.00'}\n")
        f.write(f"NUCLEAR_REP {'0.0'}\n")

        f.write(f"INITIAL_E_T {'-254.912742529660'}\n")
        f.write(f"PRINT_WALKER_INFO {'false'}\n")
        f.write(f"FORCE_SLICE  {'false'}\n")
        f.write(f"FORCE_SLICEGROUPS_BOTH {'false'}\n")
        f.write(f"FORCE_SLICEGROUPS_Y  {'false'}\n")
        f.write(f"FORCE_NGROUPS {'0'}\n")
        f.write(f"TEST_OPTSLICE {'false'}\n")
        f.write(f"MEASURE_MEM {'false'}\n")
        f.write(f"READ_WALKER_CHK {'false'}\n")
        f.write(f"WRITE_WALKER_CHK {'false'}\n")
        f.write(f"NBLK_WRITE_WALKER_CHK {'100'}\n")
        f.write(f"READ_WALKER_CHK_DIR \"{'./walker_files'}\"\n")
        f.write(f"READ_Y_DIR \"{'Y_files'}\"\n")
        f.write(f"READ_VLIST_DIR \"{'vList_files'}\"\n")
        f.write(f"READ_EXTRA_1BODY {'false'}\n")
        f.write(f"WRITE_EXTRA_1BODY {'false'}\n")
        f.write(f"FLAG_RHFBASIS {'false'}\n")
        f.write(f"FLAG_RHFTRIAL {'false'}\n")
        f.write(f"FLAG_UHFBASIS {'false'}\n")
        f.write(f"FLAG_MEASURETRIAL {'false'}\n")
        f.write(f"N_STREAMS {'32'}\n")
        f.write(f"DONT_PRECOMPUTE_Y {'false'}\n")
        f.write(f"PRINT_DATA {'false'}\n")
        f.write(f"FLAG_SMW {'true'}\n")
        f.write(f"COMPRESS_ERIS {'false'}\n")
        f.write(f"COMPRESS_ERIS_REAL {'false'}\n")
        f.write(f"COMPRESS_ERIS_THRESHOLD {'0.001'}\n")
        f.write(f"READ_COMPRESS_ERIS {'false'}\n")
        f.write(f"WRITE_COMPRESS_ERIS {'false'}\n")
        f.write(f"MPI_SPLIT_COMPRESSION {'true'}\n")
        f.write(f"FLAG_CHOLESKY_SYMM {'false'}\n")
        f.write(f"FLAG_CHOLESKY_SYMMREAL {'true'}\n")
        f.write(f"CORRELATED_SAMPLING {'false'}\n")
        f.write(f"TRANSFORM_LOCALIZED_MOS {'false'}\n")
        f.write(f"SMW_BATCH {'false'}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fcidump_file> <ezfio_filename>")
        sys.exit(1)
    else:
        fcidump_file = sys.argv[1]
        ezfio_filename = sys.argv[2]

        NORB, n_basis, n_electron, alpha_num, beta_num, n_det, afqmc_in_filename,filename = read_ezfio(ezfio_filename)
        eri = dump_fci(fcidump_file, NORB)
        mcd= core_mcd(eri, NORB, chmax=2000, tolcd=1e-5)
        eri = dump_fci(fcidump_file, NORB)
        write_mcd_core(filename, NORB, mcd)
        afqmc_in(n_basis, n_electron, alpha_num, beta_num, n_det, afqmc_in_filename,filename)
