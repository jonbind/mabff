#!/usr/bin/env python3
import os
from typing import Dict, List

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser


"""
Coarse-grained antibody topology builder.

Converts an all-atom PDB structure into a one-bead-per-residue (1bAA)
coarse-grained representation and writes a GROMACS topology (.itp) with an
elastic network model, as described in Binder et al., "Predicting
High-Concentration Antibody Dynamics: Development and Validation of a
Sequence-Specific Coarse-Grained Model".

Usage (CLI):
    python build_cg_topology.py <input.pdb> <molecule_name> <output.itp>

Workflow:
1) Coarse-graining:
   The input PDB is reduced to one 'BB' bead per residue, placed at the
   residue's geometric center, and written to <input>_cg.pdb.
2) Glycan placement:
   Not implemented; glycans are omitted in the current model iteration.
   Beads whose names start with 'G' are reserved for future glycan support
   and are explicitly skipped when building bonds and angles.
3) Topology generation:
   The coarse-grained PDB is read back and written as a GROMACS .itp file.

The model parameters (bond/angle force constants, elastic-network cutoffs and
charge assignment) are defined as named constants below for transparency.

This module exports:
- calculate_geometric_center(residue)
- list_residues(pdb_file, temp_file)        # writes coarse-grained PDB (fixed-width lines)
- preparemolecule(temp_file)                # reads the coarse-grained PDB into a DataFrame
- writemolecule(mAb_df, itp_output_file_path, molecule_name)
- run_workflow(input_pdb, molecule_name, output_itp)
"""


# -------------------------------------------------------------------
# Coarse-grained model parameters
# (see Binder et al., Materials and Methods — Coarse-Graining Scheme)
# -------------------------------------------------------------------
BEAD_NAME = "BB"                  # one backbone bead per residue

# Backbone bonds between consecutive residues in the same chain
BACKBONE_BOND_FUNCT = 1
BACKBONE_BOND_LENGTH = 0.4        # nm
BACKBONE_BOND_FORCE = 6000        # kJ mol^-1 nm^-2

# Elastic network model ("rubber bands") between residue pairs whose initial
# separation falls within [ENM_LOWER, ENM_UPPER].
ENM_FUNCT = 1
ENM_LOWER = 0.4                   # nm, lower distance cutoff
ENM_UPPER = 0.6                   # nm, upper distance cutoff
ENM_FORCE = 600.0                 # kJ mol^-1 nm^-2

# Backbone angles between consecutive residue triplets in the same chain
ANGLE_FUNCT = 10
ANGLE_THETA = 127                 # degrees
ANGLE_FORCE = 20                  # kJ mol^-1 rad^-2

# Position restraints applied under -DPOSRES (isotropic, x/y/z)
POSRES_FORCE = 1000               # kJ mol^-1 nm^-2

# Residue charge assignment. The protonation state is encoded in the residue
# name upstream (see rename_protonation_states.py); all other residues,
# including neutral protonation variants (e.g. GLH, ASH), carry zero charge.
POSITIVE_RESIDUES = {"ARG", "LYS", "HIP", "HSP"}
NEGATIVE_RESIDUES = {"ASP", "GLU"}


# -------------------------
# Coarsegrainer (Step 1)
# -------------------------

def calculate_geometric_center(residue) -> np.ndarray:
    """
    Calculate the geometric center of a residue based on the coordinates of its atoms.
    Returns a 3D numpy array [x, y, z].
    """
    coords = [atom.get_coord() for atom in residue.get_atoms()]
    if not coords:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    return np.mean(coords, axis=0)


def _format_pdb_atom_line(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    occupancy: float = 1.00,
    bfactor: float = 0.00,
    element: str = "C",
    charge: str = "",
) -> str:
    """
    Format a single ATOM line in fixed-width PDB columns.

    Columns (1-based):
    1-6   Record name  "ATOM  "
    7-11  serial       (right aligned)
    13-16 atom name    (4 chars)
    17    altLoc
    18-20 resName
    22    chainID
    23-26 resSeq
    27    iCode
    31-38 x (8.3)
    39-46 y (8.3)
    47-54 z (8.3)
    55-60 occupancy (6.2)
    61-66 tempFactor (6.2)
    77-78 element (right aligned)
    79-80 charge
    """
    record = "ATOM"
    altloc = ""
    icode = ""
    name = (name[:4]).ljust(4)  # left-justified atom name
    resname = (resname[:3]).upper()
    chain = (chain or "A")[:1]
    element = (element or "").upper()[:2]
    charge = (charge or "")[:2]

    return (
        f"{record:<6}{serial:>5d} {name}{altloc:1}"
        f"{resname:>3} {chain:1}{resseq:>4d}{icode:1}"
        f"   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{bfactor:>6.2f}"
        f"          {element:>2}{charge:>2}\n"
    )


def list_residues(pdb_file: str, temp_file: str) -> int:
    """
    Coarse-grain writer:
    For each residue, write a single 'BB' bead at the geometric center in proper PDB fixed-width format.
    Returns the number of residues written.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_file)

    count = 0
    serial = 1
    with open(temp_file, "w") as fh:
        for chain in structure.get_chains():
            chain_id = chain.get_id().strip() or "A"
            for residue in chain.get_residues():
                hetflag, resseq, _icode = residue.get_id()

                coords = [atom.get_coord() for atom in residue.get_atoms()]
                if not coords:
                    continue
                center = np.mean(coords, axis=0)
                if np.isnan(center).any():
                    continue

                resname = residue.get_resname()
                line = _format_pdb_atom_line(
                    serial=serial,
                    name=BEAD_NAME,
                    resname=resname,
                    chain=chain_id,
                    resseq=int(resseq),
                    x=float(center[0]),
                    y=float(center[1]),
                    z=float(center[2]),
                    occupancy=1.00,
                    bfactor=0.00,
                    element="C",
                    charge=""
                )
                fh.write(line)
                serial += 1
                count += 1

        fh.write("END\n")

    return count


# -------------------------
# Molprep (Step 3)
# -------------------------

def preparemolecule(temp_file: str) -> pd.DataFrame:
    """
    Read the coarse-grained PDB (fixed-width) and build the DataFrame used by writemolecule.
    Parses fixed columns to avoid whitespace issues.

    Produces columns: Atomnumber, BeadID, Resname, Chain, ResSeq, X, Y, Z, Charge
    """
    records: List[Dict] = []
    with open(temp_file, "r") as fh:
        for line in fh:
            if line.startswith("END"):
                break
            if not line.startswith("ATOM"):
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16].strip()
                resname = line[17:20].strip().upper()
                chain = line[21:22].strip() or "A"
                resseq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                # skip malformed lines
                continue

            records.append(
                {
                    "Atomnumber": serial,
                    "BeadID": name or "BB",
                    "Resname": resname,
                    "Chain": chain,
                    "ResSeq": resseq,
                    "X": x,
                    "Y": y,
                    "Z": z,
                }
            )

    if not records:
        raise ValueError(f"No ATOM records parsed from {temp_file}")

    mAb = pd.DataFrame.from_records(records)
    mAb = mAb.sort_values(["Chain", "Atomnumber"]).reset_index(drop=True)

    # Charges by residue type (protonation encoded in residue name upstream)
    mAb["Charge"] = 0.0
    mAb.loc[mAb["Resname"].isin(POSITIVE_RESIDUES), "Charge"] = 1.0
    mAb.loc[mAb["Resname"].isin(NEGATIVE_RESIDUES), "Charge"] = -1.0

    # Ensure dtypes
    mAb["Atomnumber"] = mAb["Atomnumber"].astype(int)
    mAb["BeadID"] = mAb["BeadID"].astype(str)
    mAb["Resname"] = mAb["Resname"].astype(str)
    mAb["Chain"] = mAb["Chain"].astype(str)
    for c in ("X", "Y", "Z", "Charge"):
        mAb[c] = mAb[c].astype(float)

    return mAb


def writemolecule(
    mAb: pd.DataFrame,
    itp_output_file_path: str,
    molecule_name: str,
    export_bonds_csv: str = None,
) -> None:
    """
    Writes the GROMACS ITP file from the prepared molecule DataFrame.
    Sections: [moleculetype], [atoms], [ position_restraints ], [ bonds ], Rubber bonds, [ angles ]

    If export_bonds_csv is given, the list of elastic-network bond pairs is
    additionally written to that CSV path (disabled by default).
    """
    if os.path.exists(itp_output_file_path):
        os.remove(itp_output_file_path)

    with open(itp_output_file_path, "w") as itp:
        # Header and moleculetype
        itp.write("; itp file\n\n")
        itp.write("; Coarse-grained antibody model topology\n")
        itp.write("[ moleculetype ]\n")
        itp.write("; molname    nrexcl\n")
        itp.write(f"{molecule_name}  1\n\n")

        # Atoms
        itp.write("[ atoms ]\n")
        itp.write("; id  Type  Resnr  Resn  BeadID  Cgnr  charge\n")
        mAb2 = mAb.dropna().copy()
        mAb2.loc[:, "Atomnumber"] = mAb2["Atomnumber"].astype(int)
        for _, row in mAb2.iterrows():
            atom_id = int(row["Atomnumber"])
            resname = str(row["Resname"])
            bead_id = str(row["BeadID"])
            resnr = int(row.get("ResSeq", atom_id))
            cgnr = atom_id
            charge = float(row["Charge"])
            itp.write(
                f"{atom_id:>4d}{resname:>6}{resnr:>6d}{resname:>5}{bead_id:>6}{cgnr:>6d}{charge:>8.3f}\n"
            )

        # Position restraints
        itp.write("\n[ position_restraints ]\n")
        itp.write("#ifdef POSRES\n")
        for atom_id in mAb2["Atomnumber"].astype(int).tolist():
            itp.write(
                f"{atom_id:>6d}{1:>4d}"
                f"{POSRES_FORCE:>8d}{POSRES_FORCE:>8d}{POSRES_FORCE:>8d}\n"
            )
        itp.write("#endif\n")

        # Backbone bonds (consecutive residues in same chain; skip glycan beads starting with 'G')
        itp.write("\n[ bonds ]\n")
        itp.write("; Backbone bonds\n")
        mAb_sorted = mAb2.sort_values(["Chain", "Atomnumber"]).reset_index(drop=True)
        for i in range(len(mAb_sorted) - 1):
            chain_i = str(mAb_sorted.at[i, "Chain"])
            chain_j = str(mAb_sorted.at[i + 1, "Chain"])
            bead_next = str(mAb_sorted.at[i + 1, "BeadID"])
            if (chain_i == chain_j) and (not bead_next.startswith("G")):
                a1 = int(mAb_sorted.at[i, "Atomnumber"])
                a2 = int(mAb_sorted.at[i + 1, "Atomnumber"])
                itp.write(
                    f"{a1:>6d}{a2:>6d}{BACKBONE_BOND_FUNCT:>4d}"
                    f"{BACKBONE_BOND_LENGTH:>8.3f}{BACKBONE_BOND_FORCE:>8d}\n"
                )

        # Rubber bonds between pairs 0.4-0.6 nm (excluding beads G0/G1)
        itp.write("; Rubber\n")
        mask = ~mAb2["BeadID"].isin(["G0", "G1"])
        mf = mAb2.loc[mask].copy()
        if not mf.empty and len(mf) > 1:
            coords = mf[["X", "Y", "Z"]].to_numpy(dtype=float)  # Angstrom
            diff = coords[:, None, :] - coords[None, :, :]
            dist = np.sqrt(np.sum(diff * diff, axis=-1)) / 10.0  # nm
            n = dist.shape[0]
            iu, ju = np.triu_indices(n, k=1)
            pair_dist = dist[iu, ju]
            sel = (pair_dist > ENM_LOWER) & (pair_dist < ENM_UPPER)
            iu_sel = iu[sel]
            ju_sel = ju[sel]
            d_sel = pair_dist[sel]

            if iu_sel.size > 0:
                atom_ids = mf["Atomnumber"].to_numpy(dtype=int)
                lines = []
                bonds_rows = []
                for a_idx, b_idx, d in zip(iu_sel, ju_sel, d_sel):
                    a = int(atom_ids[a_idx])
                    b = int(atom_ids[b_idx])
                    lines.append(
                        f"{a:>6d}{b:>6d}{ENM_FUNCT:>4d}{d:>8.4f}{ENM_FORCE:>8.1f}\n"
                    )
                    bonds_rows.append((a, b))
                itp.writelines(lines)
                # Optionally export the elastic-network bond pairs to CSV
                if export_bonds_csv:
                    pd.DataFrame(bonds_rows, columns=["Atom1", "Atom2"]).to_csv(
                        export_bonds_csv, index=False, header=False
                    )
            else:
                itp.write(f"; No pairs within {ENM_LOWER}-{ENM_UPPER} nm.\n")
        else:
            itp.write("; Rubber: not enough eligible beads.\n")

        # Backbone angles (triplets in same chain; skip glycan beads)
        itp.write("\n[ angles ]\n")
        itp.write("; Backbone angles\n")
        n = len(mAb_sorted)
        for i in range(n - 2):
            chain_i = str(mAb_sorted.at[i, "Chain"])
            chain_j = str(mAb_sorted.at[i + 1, "Chain"])
            chain_k = str(mAb_sorted.at[i + 2, "Chain"])
            b1 = str(mAb_sorted.at[i + 1, "BeadID"])
            b2 = str(mAb_sorted.at[i + 2, "BeadID"])
            if chain_i == chain_j == chain_k and (not b1.startswith("G")) and (not b2.startswith("G")):
                a1 = int(mAb_sorted.at[i, "Atomnumber"])
                a2 = int(mAb_sorted.at[i + 1, "Atomnumber"])
                a3 = int(mAb_sorted.at[i + 2, "Atomnumber"])
                itp.write(
                    f"{a1:>6d}{a2:>6d}{a3:>6d}"
                    f"{ANGLE_FUNCT:>6d}{ANGLE_THETA:>8d}{ANGLE_FORCE:>8d}\n"
                )


# -------------------------
# Workflow entry point
# -------------------------

def run_workflow(input_pdb: str, molecule_name: str, output_itp: str) -> None:
    """
    Workflow:
    1) Coarsegrainer: input PDB -> [base]_cg.pdb
    2) Glycfinder: (placeholder; not implemented here)
    3) Molprep: read [base]_cg.pdb -> write ITP
    """
    if not os.path.isfile(input_pdb):
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")
    base, ext = os.path.splitext(input_pdb)
    temp_cg = f"{base}_cg{ext}"

    # Step 1: coarse-grain PDB to a single BB bead per residue
    written = list_residues(input_pdb, temp_cg)
    if written == 0:
        raise RuntimeError("No residues were written during coarse-graining.")

    # Step 2: Glycan placeholder (if sugars are added, use that output instead of temp_cg)

    # Step 3: Prepare molecule and write ITP
    mAb = preparemolecule(temp_cg)
    writemolecule(mAb, output_itp, molecule_name)


def main(argv: List[str] = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Convert an all-atom antibody PDB into a one-bead-per-residue "
            "coarse-grained GROMACS topology (.itp) with an elastic network model."
        )
    )
    parser.add_argument("input_pdb", help="Input all-atom PDB file")
    parser.add_argument("molecule_name", help="Molecule name written to [ moleculetype ]")
    parser.add_argument("output_itp", help="Output GROMACS topology (.itp) file")
    args = parser.parse_args(argv)

    run_workflow(args.input_pdb, args.molecule_name, args.output_itp)
    print(f"Wrote coarse-grained topology to {args.output_itp}")


if __name__ == "__main__":
    main()
