#!/usr/bin/env python3

import numpy as np

coeffs = []  

def occ_bin(line):
    # Convert '+ - to '10'
    return ''.join('1' if s == '+' else '0' for s in line.strip().strip('|') if s in '+-')


def get_coeffs(coeff_file):
    with open(coeff_file, 'r') as t:
        for line in t:
            stripped = line.strip()
            if stripped:
                coeffs.append(float(stripped))
    print(f"Read {len(coeffs)} coefficients from text file.")


def process_occupancy(wf, output='CI_coeff.dat'):
    # Active orbital indices (right to left, 0-based)
    Act_mos = [39, 31, 30, 28, 27, 26, 25, 24]
    #[15, 14, 13, 12]  # LSB is orbital 12

    strings = []

    with open(wf, 'r') as f:
        lines = f.readlines()


    occ_lines = [
        line.strip()
        for line in lines
        if line.strip().startswith('|') and line.strip().endswith('|') and '+' in line and '-' in line
    ]

    if len(occ_lines) % 2 != 0:
        raise ValueError("Mos can only be even).")

    n_dets = len(occ_lines) // 2

    if len(coeffs) != n_dets:
        raise ValueError(
            f"Mismatch: found {n_dets} alpha-beta pairs in file, but {len(coeffs)} coefficients were given."
        )

    for i in range(0, len(occ_lines), 2):
        alpha_line = occ_lines[i]
        beta_line  = occ_lines[i + 1]

        alpha_bits = occ_bin(alpha_line)
        beta_bits  = occ_bin(beta_line)

        alpha_active = ''.join(alpha_bits[j] for j in Act_mos)
        beta_active  = ''.join(beta_bits[j] for j in Act_mos)

        alpha_bin = f"0b{alpha_active.lstrip('0') or '0'}"
        beta_bin  = f"0b{beta_active.lstrip('0') or '0'}"

        strings.append((alpha_bin, beta_bin))

    # Compose output lines
    output_lines = []
    for idx, (occa, occb) in enumerate(strings):
        coeff = coeffs[idx]
        output_lines.append(f"{occa:>20} {occb:>20} {coeff:>20.16f}")

    # Output to file or terminal
    if output:
        with open(output, 'w') as out:
            out.write('\n'.join(output_lines))
    else:
        for line in output_lines:
            print(line)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Match CI coefficients from Trexio to alpha/beta determinants read from a .wf file.")
    parser.add_argument("input", help="wavefunction file(.wf)")
    parser.add_argument("coeff_file", help="Coefficient file from trexio (.txt)")
    parser.add_argument("-o", "--output", help="path/name to the output file you want to write", default="CI_coeff.dat")
    args = parser.parse_args()

    get_coeffs(args.coeff_file)
    process_occupancy(args.input, output=args.output)


if __name__ == "__main__":
    main()

