"""
Rename protonation-state residue names from a PDB2PQR/PROPKA-prepared structure
to the standard names expected by the coarse-graining script (mAb_prepv3.py).

Mappings applied:
    HI+  -> HIP   (protonated, positively charged histidine)
    CSS  -> CYS   (disulfide-bonded cysteine)
    GL0  -> GLH   (neutral, protonated glutamate)
    AS0  -> ASH   (neutral, protonated aspartate)

Residue identity is tracked by (chain, resSeq, iCode) so that all atoms of a
given residue are renamed consistently.

Usage:
    python rename_charged_species_in_pdb.py <input.pdb> <output.pdb>
"""

import sys


def rename_pdb_residues(input_pdb_file, output_pdb_file):
    # Map residues to new 3-letter names using a unique residue identifier
    # that includes chain, residue sequence number, and insertion code.
    rename_dict = {}

    def get_res_id(line):
        # PDB columns (1-based):
        # chainID: col 22  -> line[21]
        # resSeq:  cols 23-26 -> line[22:26]
        # iCode:   col 27  -> line[26]
        chain = line[21]
        res_seq = line[22:26].strip()
        i_code = line[26]  # may be ' ' or a letter
        return (chain, res_seq, i_code)

    # First pass: collect residues that need renaming based on existing names
    with open(input_pdb_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                res_id = get_res_id(line)

                if res_name == 'HI+':
                    # Rename protonated histidine to HIP
                    rename_dict[res_id] = 'HIP'

    # Second pass: write out with renames applied
    with open(input_pdb_file, 'r') as f, open(output_pdb_file, 'w') as out_f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                res_id = get_res_id(line)

                if res_id in rename_dict:
                    new_name = rename_dict[res_id]
                    # Overwrite residue name (cols 18-20, 0-based 17:20), keep formatting
                    line = line[:17] + f"{new_name:>3}" + line[20:]
                elif res_name == "CSS":
                    line = line[:17] + "CYS" + line[20:]
                elif res_name == "GL0":
                    line = line[:17] + "GLH" + line[20:]
                elif res_name == "AS0":
                    line = line[:17] + "ASH" + line[20:]

            out_f.write(line)


def main():
    if len(sys.argv) != 3:
        print("Usage: python rename_charged_species_in_pdb.py input.pdb output.pdb")
        sys.exit(1)
    input_pdb_file = sys.argv[1]
    output_pdb_file = sys.argv[2]
    rename_pdb_residues(input_pdb_file, output_pdb_file)


if __name__ == "__main__":
    main()