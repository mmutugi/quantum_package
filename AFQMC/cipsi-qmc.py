#!/usr/bin/env python3

try:
    import trexio
    import numpy as np
    from scipy.io import FortranFile
    import sys
    import time

except:
    print("TREXIO is not installed. Please install it first.")
    raise ImportError

import numpy

def get_afqmc_data(trexio_filename, verbose: bool = True) -> None:
    trexio_file = trexio.File(trexio_filename, "r", trexio.TREXIO_AUTO)

    mo_num = trexio.read_mo_num(trexio_file)
    if verbose:
        print(f"mo_num = {mo_num} # MOs")

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
        print(f"chol_num = {chol_num} # Cholesky vectors")

    chol = numpy.zeros((mo_num * mo_num, chol_num))
    BUFFER_SIZE = 1000000
    offset = 0
    eof = False
    while not eof:
        indices, values, nread, eof = trexio.read_mo_2e_int_eri_cholesky(
            trexio_file, offset, BUFFER_SIZE
        )
        offset += nread
        for l in range(nread):
            i, j, k = indices[l]
            integral = values[l]
            chol[i * mo_num + j, k] = integral
    L = chol
    if verbose:
        print(f"Read Cholesky vectors")

    #%_______________________________________________________
    #BUFFER_SIZE = 1000000
    #offset = 0
    #eri = np.zeros([mo_num, mo_num, mo_num, mo_num])
    #eof = False
    #fcidump_threshold = 1e-10
    #integral_eof = False
    #while not integral_eof:
     #   indices, values, nread, integral_eof = trexio.read_mo_2e_int_eri(
      #      trexio_file, offset, BUFFER_SIZE
       ##offset += nread
        #for integral in range(nread):
         #   val = values[integral]
          #  if np.abs(val) < fcidump_threshold:
           #     continue
            #i, j, k, l = indices[integral]

            #eri[i,j,k,l] = val
            #eri[k,l,i,j] = val
            #eri[j,i,l,k] = val
            #eri[l,k,j,i] = val
            #eri[j,i,k,l] = val
            #eri[l,k,i,j] = val
            #eri[i,j,l,k] = val
            #eri[k,l,j,i] = val
            #chem_i, chem_j, chem_k, chem_l = i, j, k, l

            #for p, q in ((i, j), (j, i)):
             #   for r, s in ((k, l), (l, k)):
              #      eri[p, q, r, s] = val
               #     eri[r, s, p, q] = val
                #    print(f"ERI: {p} {q} {r} {s} {val}")
    

    ##fcidumpfile


    determinants = trexio.read_determinant_list(trexio_file, 0, ndet)[0]
    nint = trexio.get_int64_num(trexio_file)

    binstr = []
    occa = []
    occb = []

    start = time.time()

    for d in determinants:
        alpha_bits = d[:nint]
        beta_bits = d[nint:]

        occa.append(trexio.to_orbital_list(nint, alpha_bits)) 
        occb.append(trexio.to_orbital_list(nint, beta_bits))

        print (f"time taken to read determinant {len(binstr)}: {time.time() - start:.2f} seconds")

        # Convert each 64-bit integer to binary and concatenate

        alpha_int = sum(x << (64 * i) for i, x in enumerate(alpha_bits[::-1]))
        beta_int = sum(x << (64 * i) for i, x in enumerate(beta_bits[::-1]))

        print (f"time for the conversion of determinant {len(binstr)}: {time.time() - start:.2f} seconds")

        alpha_bin = bin(alpha_int)
        beta_bin  = bin(beta_int)

        binstr.append((alpha_bin, beta_bin))

    occa = numpy.array(occa)
    occb = numpy.array(occb)
    binstr = numpy.array(binstr, dtype=object)

    start_write = time.time()

    # Remove consistently occupied orbitals from the right
    n_frozen_orbs = 0
    while True:
        try:
        # Reverse strings and slice from n_frozen_orbs (to check last orbitals)
            if ('0' in [i[0][::-1][n_frozen_orbs] for i in binstr]) or ('0' in [i[1][::-1][n_frozen_orbs] for i in binstr]):
                break
            else:
                n_frozen_orbs += 1
        except IndexError:
            break

    print(f"Found {n_frozen_orbs} always occupied orbitals – will drop from active space")
    act_orbs = mo_num - n_frozen_orbs
    act_electrons = nup + ndn - 2 * n_frozen_orbs
    act_alpha = nup - n_frozen_orbs + 1
    beta_act = ndn - n_frozen_orbs - 1

# Apply the truncation to the bitstrings
    binstr = [(a[:-n_frozen_orbs] if n_frozen_orbs > 0 else a,
           b[:-n_frozen_orbs] if n_frozen_orbs > 0 else b) for a, b in binstr]


    with open ('CI_coeff.dat','w') as f:
        for i in range(len(ci_coeffs)):
            coeff = ci_coeffs[i]
            alpha_bin, beta_bin = binstr[i]
            f.write(f"{alpha_bin:<30} {beta_bin :<30} {-coeff: .10f}\n")

            print(f"Writing CI coefficient {i+1}/{len(ci_coeffs)}: {time.time() - start_write:.2f} seconds")


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

                #alpha_orb = trexio.to_orbital_list(nint, determinants[i][:nint])
                #alpha_orb = [x+1 for x in alpha_orb]
                #alpha_orb_str = ' '.join([str(x) for x in alpha_orb])
                #beta_orb = trexio.to_orbital_list(nint, determinants[i][nint:])
                #beta_orb = [x+1 for x in beta_orb]
                #beta_orb_str = ' '.join([str(x) for x in beta_orb])
                #coeff = ci_coeffs[i]
                #cas_file.write(f"{alpha_orb_str} {beta_orb_str} {coeff:16f}\n")
            cas_file.write('\n')
            cas_file.write('multidet_ampl\n')
            cas_file.write('#             amplitude  # ndets   det_tot_weight\n')


            for i in range(ndet):
                cas_file.write(f"{ci_coeffs[i]:20.16f} 0.0  #{i:20}\t{wt_list[i]:20.16f}\n")
            cas_file.write('\n')
            cas_file.write('multidet_type 1\n')
            cas_file.write(f'npsitdet {ndet}\n')
            cas_file.write('# --- END TRIAL WAVE FUNCTION --- simple cut\n')

    with open('MCSCF_MOs.dat','w') as mcsf_file:
        s1 = np.identity(mo_num)
        for x in range(mo_num):
            for y in range(mo_num):
                mcsf_file.write(f"{s1[x][y]: 20.16f}")
            mcsf_file.write('\n')

    with open('ROHF.dat','w') as ROHF_file:
        identity = np.identity(mo_num)
        ROHF_file.write('%9s %9s #ROHF orbitals \n' % (mo_num,mo_num))
        ROHF_file.write('\n')
        for ip in range(mo_num):
            for jp in range(mo_num):
                ROHF_file.write('%20.16f \n' % (identity[ip,jp]))
            ROHF_file.write('\n')
        ROHF_file.write('%9s %9s #ROHF orbitals \n' % (mo_num,mo_num))
        ROHF_file.write('\n')
        for ip in range(mo_num):
            for jp in range(mo_num):
                ROHF_file.write('%20.16f \n' % (identity[ip,jp]))
            ROHF_file.write('\n')

    with open('one_body_gms','w') as one_body_file:
        nbasis = hcore.shape[0]
        identity = np.identity(nbasis)
        one_body_file.write('%9s %9s # \n' % (nbasis,nbasis*nbasis))

        one_body_file.write('Overlap \n')
        for ip in range(nbasis):
            for jp in range(nbasis):
                one_body_file.write(' %9s %9s % .12f \n' % (ip+1,jp+1,s1[ip,jp]))

        one_body_file.write('\nCore Hamiltonian \n')
        for ip in range(nbasis):
            for jp in range(nbasis):
                one_body_file.write(' %9s %9s % .12f \n' % (ip+1,jp+1,hcore[ip,jp]))
            #one_body_file.write('\n')

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
        "N_ACT_ORB": act_orbs,
        "N_ACT_ELEC": act_electrons,
        "N_ACT_UP": act_alpha,
        "N_ACT_DN": beta_act,



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


def afqmc_in(NORB,nup, N_ACT_ORB, N_ACT_UP, N_ACT_DN,  ndn, ndet, afqmc_in_filename = 'afqmc.in'):
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
        f.write(f"N_ACT_ORBS {N_ACT_ORB}\n")
        f.write(f"N_ACT_UP {N_ACT_UP}\n")
        f.write(f"N_ACT_DN {N_ACT_DN}\n")
        f.write(f"N_DET {ndet}\n")

        f.write(f"N_WLK {'1000'}\n")
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


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <trexio_file> ")
        sys.exit(1)

    get_afqmc_data(sys.argv[1])
    L = get_afqmc_data(sys.argv[1])['L']
    NORB = get_afqmc_data(sys.argv[1])['NORB']
    #mcd = core_mcd(h2e, NORB, chmax=2000, tolcd=1e-5)
    write_mcd_core(NORB, L)
    nup, ndn, ndet = get_afqmc_data(sys.argv[1])['nup'], get_afqmc_data(sys.argv[1])['ndn'], get_afqmc_data(sys.argv[1])['ndet']
    afqmc_in(NORB,nup, ndn, ndet, afqmc_in_filename = 'afqmc.in')

