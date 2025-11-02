#!/usr/bin/env python3
"""
AFQMC Input Generator from TREXIO Output
========================================

This script extracts and converts quantum chemistry data from a TREXIO-formatted directory
(ASCII or HDF5) into the required input formats for RICE afqmc. Generates:

- ROHF.dat            : Molecular orbital (MO) coefficient matrix.
- MCSCF_MOs.dat       : Same MO matrix, used for MCSCF-type trial wavefunctions.
- one_body_gms.dat    : One-body integrals (core Hamiltonian and AO overlap).
- CI_coeff.dat        : CI coefficients for multi-determinant trial wavefunctions.
- CAS_WF_rcas         : Multi-determinant trial wavefunction in real-space AFQMC format.
- afqmc.in            : AFQMC runtime input file.
- V2b_AO_cholesky.mat : Two-electron integrals (Cholesky-decomposed) in binary Fortran format.

Main Functional Components:
---------------------------

- `read_mo_coeff`:
    Parses `mo.txt` from TREXIO to extract the number of MOs and the MO coefficient matrix.

- `read_overlap`, `read_hamiltonian`:
    Extracts the AO overlap and one-electron Hamiltonian integrals from `ao_1e_int.txt`.

- `write_rohf_file`, `write_mcsf_file`:
    Writes MO coefficients to standard text formats used by AFQMC codes.

- `get_afqmc_data`:
    Uses the TREXIO Python API to extract multi-determinant trial wavefunctions, CI coefficients,
    and Cholesky vectors from HDF5 TREXIO data.

- `write_one_body_gms`, `write_mcd_core`, `afqmc_in`:
    Write out AFQMC-specific formatted input files (text and Fortran binary).

Usage:
------

Run from the command line by providing the path to the TREXIO directory (e.g., `test.trexio`):

    python full.py path/to/test.trexio

By default, verbose mode is enabled. To disable verbose output:

    python full.py path/to/test.trexio --no-verbose

Assumptions:
------------
- The TREXIO directory must include:
    - `mo.txt`              : containing `mo_num` and `mo_coefficient`
    - `ao_1e_int.txt`       : containing `ao_1e_int_overlap` and `ao_1e_int_core_hamiltonian`

Setup TRexio:
--------------
qp set_file xyz.ezfio
qp set trexio backend 1
qp set trexio trexio_file xyz.trexio
qp set trexio export_basis 1
qp set trexio export_ao_one_e_ints 1
qp set trexio export_ao_two_e_ints 1
qp set trexio export_ao_two_e_ints_cholesky 1
qp set trexio export_mo_one_e_ints 1
qp set trexio export_mo_two_e_ints 1
qp set trexio export_mo_two_e_ints_cholesky 1
qp set trexio export_rdm 1
qp run export_trexio

Author:
-------
Mark Munyi

"""


try:
    import trexio
    import numpy as np
    from scipy.io import FortranFile
    import sys
    import argparse
    import os

except:
    print("TREXIO is not installed. Please install it first.")
    raise ImportError

def read_mo_coeff(mo_txt_file: str, verbose: bool = True) -> np.ndarray:
    with open(mo_txt_file, 'r') as f:
        lines = f.readlines()

    # Find mo_num
    mo_num = None
    for line in lines:
        if line.strip().startswith("mo_num "):  # Avoid mo_num_isSet
            parts = line.strip().split()
            if len(parts) == 2:
                mo_num = int(parts[1])
                if verbose:
                    print(f"Found mo_num: {mo_num}")
                break
    if mo_num is None:
        raise ValueError("mo_num not found in TREXIO file.")

    # Find start of mo_coefficient block
    start_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "mo_coefficient":
            start_idx = i + 1
            if verbose:
                print(f"Found start of mo_coefficient at line {start_idx}")
            break
    if start_idx is None:
        raise ValueError("mo_coefficient block not found.")

    # Read mo_num * mo_num coefficients
    num_coeffs = mo_num * mo_num
    coeffs = []
    for line in lines[start_idx:]:
        if line.strip() == "" or not any(c.isdigit() for c in line):
            continue  # skip empty or non-numeric lines
        try:
            val = float(line.strip())
            coeffs.append(val)
            if len(coeffs) >= num_coeffs:
                break
        except ValueError:
            continue

    if len(coeffs) != num_coeffs:
        raise ValueError(f"Expected {num_coeffs} MO coefficients, found {len(coeffs)}.")

    coeff_matrix = np.array(coeffs).reshape((mo_num, mo_num))
    return coeff_matrix


def write_rohf_file(filename: str, mo_coeff: np.ndarray) -> None:
    mo_num = mo_coeff.shape[0]
    identity = np.identity(mo_num)

    with open(filename, 'w') as f:
        f.write(f"{mo_num:9d} {mo_num:9d} # ROHF orbitals\n\n")
        for row in mo_coeff:
            for val in row:
                f.write(f"{val:20.16f}\n")
            f.write("\n")
        f.write(f"{mo_num:9d} {mo_num:9d} # ROHF orbitals\n\n")
        for row in mo_coeff:
            for val in row:
                f.write(f"{val:20.16f}\n")
            f.write("\n")
    return mo_num


def write_mcsf_file(filename: str, mo_coeff: np.ndarray) -> None:
    with open(filename, 'w') as f:
        for row in mo_coeff:
            for val in row:
                f.write(f"{val:20.16f}")
            f.write("\n")

def read_overlap(hamiltonian_file: str,mo_num, verbose: bool = True) -> np.ndarray:
    with open(hamiltonian_file, 'r') as f:
        lines = f.readlines()
    ##find the overlap information
    start_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "ao_1e_int_overlap":
            start_idx = i + 1
            if verbose:
                print(f"Found start the overlap info line {start_idx}")
            break
    if start_idx is None:
        raise ValueError("Cannot find overlap info.")
    
    num_coeffs = mo_num * mo_num
    overlap = []
    for line in lines[start_idx:]:
        if line.strip() == "" or not any(c.isdigit() for c in line):
            continue  # skip empty or non-numeric lines
        try:
            val = float(line.strip())
            overlap.append(val)
            if len(overlap) >= num_coeffs:
                break
        except ValueError:
            continue
    if len(overlap) != num_coeffs:
        raise ValueError(f"Expected {num_coeffs} overlaps, found {len(overlap)}.")

    overlap_matrix = np.array(overlap).reshape((mo_num, mo_num))
    return overlap_matrix

def read_hamiltonian(hamiltonian_file: str,mo_num, verbose: bool = True) -> np.ndarray:

    with open(hamiltonian_file, 'r') as f:
        lines = f.readlines()
    ##find the core hamiltonian information

    start_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "ao_1e_int_core_hamiltonian":
            start_idx = i + 1
            if verbose:
                print(f"Found start core hamiltonian {start_idx}")
            break
    if start_idx is None:
        raise ValueError("no core hamiltonian.")
    
    num_coeffs = mo_num * mo_num
    core_hamiltonian = []
    for line in lines[start_idx:]:
        if line.strip() == "" or not any(c.isdigit() for c in line):
            continue  # skip empty or non-numeric lines
        try:
            val = float(line.strip())
            core_hamiltonian.append(val)
            if len(core_hamiltonian) >= num_coeffs:
                break
        except ValueError:
            continue
    if len(core_hamiltonian) != num_coeffs:
        raise ValueError(f"Expected {num_coeffs} ents, found {len(core_hamiltonian)}.")

    hamiltonian_matrix = np.array(core_hamiltonian).reshape((mo_num, mo_num))
    return hamiltonian_matrix

def write_one_body_gms(filename: str, overlap_matrix, hamiltonian_matrix, mo_num: int) -> None:
    with open(filename, 'w') as f:
        f.write('%9s %9s # \n' % (mo_num,mo_num*mo_num))
        f.write('Overlap \n')
        for ip in range(mo_num):
            for jp in range(mo_num):
                f.write(' %9s %9s % .12f \n' % (ip+1,jp+1,overlap_matrix[ip,jp]))
        f.write('\nCore Hamiltonian \n')
        for ip in range(mo_num):
            for jp in range(mo_num):
                f.write(' %9s %9s % .12f \n' % (ip+1,jp+1,hamiltonian_matrix[ip,jp]))
    return mo_num


def write_orbitals(trexio_dir: str, verbose: bool = True) -> None:

    mo_txt_file = os.path.join(trexio_dir, "mo.txt")
    if not os.path.isfile(mo_txt_file):
        raise FileNotFoundError(f"Expected mo.txt in {trexio_dir}, but not found.")
    mo_coeff = read_mo_coeff(mo_txt_file, verbose)
    if verbose:
        print(f"MO coefficient matrix shape: {mo_coeff.shape}")

    write_rohf_file("ROHF.dat", mo_coeff)
    write_mcsf_file("MCSCF_MOs.dat", mo_coeff)

    if verbose:
        print("Files written: ROHF.dat, MCSCF_MOs.dat")

def get_afqmc_data(trexio_dir, mo_num, verbose: bool = True) -> None:
    trexio_file = trexio.File(trexio_dir, "r", trexio.TREXIO_AUTO)


    ndet = trexio.read_determinant_num(trexio_file)
    if verbose:
        print(f"{ndet} determinants")

    nup = trexio.read_electron_up_num(trexio_file)
    ndn = trexio.read_electron_dn_num(trexio_file)
    if verbose:
        print(f"{nup} up- and {ndn} down-spin electrons")

    hcore = trexio.read_mo_1e_int_core_hamiltonian(trexio_file)
    if verbose:
        print("Core Hamiltonian found")

    ci_coeffs = trexio.read_determinant_coefficient(trexio_file, 0, ndet)[0]
    if verbose:
        print(f"Read CI coefficients")

    e0 = trexio.read_nucleus_repulsion(trexio_file)
    if verbose:
        print(f"Nucleus repulsion energy: {e0}")


    chol_num = trexio.read_mo_2e_int_eri_cholesky_num(trexio_file)
    if verbose:
        print(f"chol_num = {chol_num}") 
    chol = np.zeros((mo_num, mo_num, chol_num))

    BUFFER_SIZE = 1000000
    offset = 0
    eof = False
    while not eof:
        indices, values, nread, eof = trexio.read_mo_2e_int_eri_cholesky(
            trexio_file, offset, BUFFER_SIZE
        )
        offset += nread
        for l, integral in enumerate(values):
            i, j, k = indices[l]
            chol[i, j, k] = integral
    L = chol.reshape(mo_num * mo_num, chol_num)
    if verbose:
        print(f"Read Cholesky vectors")
        
    determinants = trexio.read_determinant_list(trexio_file, 0, ndet)[0]
    nint = trexio.get_int64_num(trexio_file)

    binstr = []
    occa = []
    occb = []

    for d in determinants:
        alpha_bits = d[:nint]
        beta_bits = d[nint:]

        occa.append(trexio.to_orbital_list(nint, alpha_bits)) 
        occb.append(trexio.to_orbital_list(nint, beta_bits))

        # Convert each 64-bit integer to binary and concatenate

        alpha_int = sum(x << (64 * i) for i, x in enumerate(alpha_bits[::-1]))
        beta_int = sum(x << (64 * i) for i, x in enumerate(beta_bits[::-1]))

        alpha_bin = bin(alpha_int)
        beta_bin  = bin(beta_int)

        binstr.append((alpha_bin, beta_bin))

    occa = np.array(occa)
    occb = np.array(occb)
    binstr = np.array(binstr, dtype=object)

    with open ('CI_coeff.dat','w') as f:
        for i in range(len(ci_coeffs)):
            coeff = ci_coeffs[i]
            alpha_bin, beta_bin = binstr[i]
            f.write(f"{alpha_bin:<30} {beta_bin :<30} {coeff: 16f}\n")

    wt_list  = np.zeros(len(determinants))
    csum = 0.0

    for i in range(len(determinants)):
        csum +=ci_coeffs[i]**2
        wt_list[i] = csum

        format_det = lambda z: '  '.join([str(i[0]-1) for i in filter(lambda x : x[1] == '1', enumerate(z))])

        with open('CAS_WF_rcas','w')as cas_file:
            cas_file.write("# --- BEGIN TRIAL WAVE FUNCTION --- simple cut\n")
            cas_file.write(f"# Total number of determinants = {ndet}\n")
            cas_file.write(f"# Total weight                 = {wt_list[ndet-1]}\n")
            cas_file.write('multidet_cfg\n')

            for i in range(ndet):
                cas_file.write(f"{format_det(binstr[i][0])}\t{format_det(binstr[i][1])}\t#\t{ci_coeffs[i]: 20.16f}\n")
            cas_file.write('\n')
            cas_file.write('multidet_ampl\n')
            cas_file.write('#             amplitude  # ndets   det_tot_weight\n')


            for i in range(ndet):
                cas_file.write(f"{ci_coeffs[i]:20.16f} 0.0  #{i:20}\t{wt_list[i]:20.16f}\n")
            cas_file.write('\n')
            cas_file.write('multidet_type 1\n')
            cas_file.write(f'npsitdet {ndet}\n')
            cas_file.write('# --- END TRIAL WAVE FUNCTION --- simple cut\n')
    if verbose:
        print(f"Read determinants")
        #print("binstr: ", binstr)
        #print("occa: ", occa)
        #print(f"First determinant: {determinants[0]}")
       # print(f"First determinant alpha (bin): {binstr[0][0]}")
        #print(f"First determinant beta (bin): {binstr[0][1]}")

    result = {
        "norb": mo_num,
        "ndet": ndet,
        "nup": nup,
        "ndn": ndn,
        "hcore": hcore,
        "ci_coeffs": ci_coeffs,
        "determinant": determinants[0],
        "chol": chol,
        "L": L[:, :chol_num],
        "NORB": mo_num,
    }
    return result



#def core_mcd(chol, NORB, chmax=2000, tolcd=1e-5):
    L = np.zeros((NORB * NORB, chmax))
    Dmax = np.zeros((NORB * NORB))
    ng = 0

    V = h2e.reshape(NORB**2, NORB**2) / 2
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
    return L[:, :ng], NORB


def write_mcd_core(NORB, L, filename ='V2b_AO_cholesky.mat'):
    with FortranFile(filename, 'w') as g:
        g.write_record(2)
        print(type(L))
        g.write_record(np.array([NORB * (NORB + 1) // 2, L.shape[1]], dtype=np.int32))
        Lmat = np.zeros((NORB, NORB))
        Llst = np.zeros((NORB * (NORB + 1) // 2))

        for gamma in range(L.shape[1]):
            Lmat = L[:, gamma].reshape((NORB, NORB))
            counter = 0
            for ib in range(NORB):
                for jb in range(ib, NORB):
                    Llst[counter] = Lmat[ib, jb]
                    counter += 1
            g.write_record(Llst[:])
        return filename


def afqmc_in(NORB,nup, ndn, ndet, afqmc_in_filename = 'afqmc.in'):
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
        f.write(f"M_BASIS {NORB}\n")
        n_electron = nup + ndn
        f.write(f"N_ELEC {n_electron}\n")
        f.write(f"ENERGY_FC {'0.0'}\n")
        f.write(f"N_UP {nup}\n")
        f.write(f"N_DN {ndn}\n")
        f.write(f"N_ACT_ORBS {NORB}\n")
        f.write(f"N_ACT_UP {nup}\n")
        f.write(f"N_ACT_DN {ndn}\n")
        f.write(f"N_DET {ndet}\n")

        f.write(f"N_WLK {'200'}\n")
        f.write(f"N_BLK {'1'}\n")
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
        f.write(f"FLAG_MEASURETRIAL {'true'}\n")
        f.write(f"N_STREAMS {'32'}\n")
        f.write(f"DONT_PRECOMPUTE_Y {'false'}\n")
        f.write(f"PRINT_DATA {'true'}\n")
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


def main():
    parser = argparse.ArgumentParser(description="Extract info for AFQMC from a trexio file.")
    parser.add_argument("trexio_filename", type=str, help="Path to TREXIO folder or file.")
    parser.add_argument("--no-verbose", dest="verbose", action="store_false", help="Disable verbose output.")
    parser.set_defaults(verbose=True)

    args = parser.parse_args()

    write_orbitals(args.trexio_filename, args.verbose)

    mo_txt_file = os.path.join(args.trexio_filename, "mo.txt")
    mo_coeff = read_mo_coeff(mo_txt_file, verbose=args.verbose)
    mo_num = mo_coeff.shape[0]

    hamiltonian_file = os.path.join(args.trexio_filename, "ao_1e_int.txt")
    overlap = read_overlap(hamiltonian_file, mo_num, args.verbose)
    hcore = read_hamiltonian(hamiltonian_file, mo_num, args.verbose)
    write_one_body_gms("one_body_gms", overlap, hcore, mo_num)

    afqmc_data = get_afqmc_data(args.trexio_filename, mo_num, args.verbose)

    afqmc_in(
        afqmc_data["NORB"],
        afqmc_data["nup"],
        afqmc_data["ndn"],
        afqmc_data["ndet"],
        afqmc_in_filename="afqmc.in"
    )
    write_mcd_core(afqmc_data["NORB"], afqmc_data["L"])

    if args.verbose:
        print("AFQMC input generation complete.")

if __name__ == "__main__":
    main()
