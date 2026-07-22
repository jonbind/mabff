# Sequence-Specific Coarse-Grained Antibody Model

Tools for converting all-atom monoclonal antibody (mAb) structures into a
one-bead-per-residue (1bAA) coarse-grained (CG) representation and generating
GROMACS topologies with an elastic network model, as described in:

> Binder J., Kalra P., Holzinger J., Thilmann C., Huelsmeyer M., Siedler M.,
> Schuetz A., Ploetz E., Lamb D. C., Friess W.
> *Predicting High-Concentration Antibody Dynamics: Development and Validation
> of a Sequence-Specific Coarse-Grained Model.* (manuscript in preparation)

If you use this code, please cite the paper above.

---

## Overview

The pipeline turns an all-atom PDB structure into a CG GROMACS topology in
three steps:

1. **Protonation assignment** (external) — assign pH-dependent protonation
   states to titratable residues, e.g. with [PDB2PQR](https://www.poissonboltzmann.org/)
   / PROPKA at the target formulation pH.
2. **Residue renaming** (`rename_protonation_states.py`) — convert the
   protonation-state residue names produced upstream to the standard names the
   coarse-graining script expects.
3. **Coarse-graining + topology** (`build_cg_topology.py`) — map each residue to a
   single `BB` bead at its geometric center and write a GROMACS `.itp` file
   containing atoms, charges, position restraints, backbone bonds, an elastic
   network, and backbone angles.

A simple graphical front-end (`gui.py`) wraps step 3.

---

## Repository contents

| File | Description |
|------|-------------|
| `scripts/rename_protonation_states.py` | Renames protonation-state residues (HI+→HIP, CSS→CYS, GL0→GLH, AS0→ASH). |
| `scripts/build_cg_topology.py` | Core coarse-graining and GROMACS topology builder (CLI + importable module). |
| `scripts/gui.py` | Optional PyQt5 GUI wrapping the coarse-graining workflow. |
| `simulation_files/` | Example GROMACS `.mdp` run parameters, force field, and solvent files. |
| `requirements.txt` | Python dependencies. |

---

## Installation

Requires Python 3.8+.

```bash
pip install -r requirements.txt
```

---

## Usage

### 1. Assign protonation states (external)

Prepare your all-atom structure with pH-dependent protonation, for example:

```bash
pdb2pqr --ff=... --with-ph=5.5 input.pdb input_ph.pdb
```

### 2. Rename charged species

```bash
python scripts/rename_protonation_states.py input_ph.pdb input_renamed.pdb
```

### 3. Build the coarse-grained topology

Command line:

```bash
python scripts/build_cg_topology.py input_renamed.pdb MOL output.itp
```

This writes:
- `input_renamed_cg.pdb` — the coarse-grained structure (one `BB` bead per residue)
- `output.itp` — the GROMACS topology

Or launch the GUI:

```bash
python scripts/gui.py
```

---

## Building a simulation system

Once you have a coarse-grained antibody topology (`mAb.itp`) and structure
(`mAb_cg.pdb`), the following [GROMACS](https://www.gromacs.org/) workflow builds
a solvated, buffered box at a target antibody concentration and runs
equilibration and production. The commands below use `gmx_mpi`; replace with
`gmx` for a serial build. The GPU flags shown in the run steps are optional.

The required inputs from this repository are `mabff213.itp` (force field),
`water.itp` / `water.gro` (solvent), and the `.mdp` run parameters in
[`simulation_files/`](simulation_files/). You additionally need `minimization.mdp`,
an ions topology (`ions.itp`), and single-bead structures for the counter-ions
and buffer species (e.g. `NA.pdb`, `CL.pdb`, `HST.pdb`).

### 1. Set up the box and add antibody copies

```bash
# center the CG antibody in a cubic box of edge L (nm)
gmx_mpi editconf -f mAb_cg.pdb -d 1.0 -o mAb.gro
gmx_mpi editconf -f mAb.gro -c -bt triclinic -box L L L -o mAb.gro

# position restraints for the antibody (used under -DPOSRES)
echo 1 | gmx_mpi genrestr -f mAb.gro -o posre.itp

# replicate to N copies to reach the target concentration
gmx_mpi insert-molecules -f mAb.gro -ci mAb.gro -o system.gro \
        -nmol <N-1> -try 5000 -selrpos atom
```

The number of antibody molecules `N` for a target concentration `c` (mg/mL) in a
cubic box of edge `L` (nm) is `N = c · L³ · 1e-21 / (MW · 1.6605e-18)`, with a
molecular weight `MW ≈ 148` kDa.

### 2. Write the topology

Create `top.top`:

```
#include "mabff213.itp"
#include "ions.itp"
#include "mAb.itp"
#ifdef POSRES
#include "posre.itp"
#endif
#include "water.itp"

[ system ]
antibody system

[ molecules ]
MOL   <N>
```

### 3. Vacuum minimization and net-charge check

```bash
gmx_mpi grompp -p top.top -f minimization.mdp -c system.gro \
        -o min1.tpr -maxwarn 2 2> charge.log
gmx_mpi mdrun -deffnm min1 -v -ntomp 8
```

`grompp` reports the system's net charge in `charge.log`
(`System has non-zero total charge: ...`).

### 4. Neutralize and add buffer/salt ions

```bash
# neutralize the net charge Q: add NA (if Q < 0) or CL (if Q > 0), |Q| copies
gmx_mpi insert-molecules -f min1.gro -ci NA.pdb -o neutralized.gro \
        -nmol <|Q|> -try 5000 -selrpos atom

# add buffer ions, n_buf of each (e.g. histidine buffer: HST+ / CL-)
gmx_mpi insert-molecules -f neutralized.gro -ci HST.pdb -o buf1.gro \
        -nmol <n_buf> -try 5000 -selrpos atom
gmx_mpi insert-molecules -f buf1.gro -ci CL.pdb -o buf2.gro \
        -nmol <n_buf> -try 5000 -selrpos atom
```

Append the added ions to the `[ molecules ]` section of `top.top`, e.g.
`NA <|Q|>`, `HST <n_buf>`, `CL <n_buf>`. The ion count for a target
concentration `C` (mol/L) is `n = C · V · N_A`, with box volume `V = L³ · 1e-24` L.

### 5. Second minimization and solvation

```bash
gmx_mpi grompp -p top.top -f minimization.mdp -c buf2.gro -o min2.tpr
gmx_mpi mdrun -deffnm min2 -v -ntomp 8

gmx_mpi solvate -cp min2.gro -cs water.gro -radius 0.44 -o solvated.gro -p top.top
```

### 6. Equilibration and production

```bash
# energy minimization
gmx_mpi grompp -p top.top -c solvated.gro -f minimization.mdp -r solvated.gro -o em.tpr
gmx_mpi mdrun -deffnm em -v

# NVT equilibration
gmx_mpi grompp -p top.top -c em.gro -f nvt.mdp -r em.gro -o nvt.tpr -maxwarn 3
gmx_mpi mdrun -deffnm nvt -v -nb gpu -bonded gpu -update gpu -pme gpu -pin on

# NPT equilibration
gmx_mpi grompp -p top.top -c nvt.gro -f npt.mdp -r nvt.gro -t nvt.cpt -o npt.tpr -maxwarn 3
gmx_mpi mdrun -deffnm npt -v -nb gpu -bonded gpu -update gpu -pme gpu -pin on

# production
gmx_mpi grompp -p top.top -c npt.gro -f dynamic.mdp -r npt.gro -t npt.cpt -o dynamic.tpr -maxwarn 3
gmx_mpi mdrun -deffnm dynamic -v -nb gpu -bonded gpu -update gpu -pme gpu -pin on
```

The production trajectory (`dynamic.xtc`) is the input for downstream dynamics
and viscosity analysis.

---

## Model parameters

The CG model parameters are defined as named constants at the top of
`build_cg_topology.py` for transparency and reproducibility:

| Parameter | Value | Description |
|-----------|-------|-------------|
| Backbone bond length | 0.4 nm | Equilibrium bond between consecutive residues |
| Backbone bond force constant | 6000 kJ·mol⁻¹·nm⁻² | |
| Elastic network cutoffs | 0.4–0.6 nm | Residue pairs restrained if their initial separation falls in this range |
| Elastic network force constant | 600 kJ·mol⁻¹·nm⁻² | Harmonic "rubber band" restraints |
| Backbone angle | 127° | Equilibrium angle for residue triplets |
| Backbone angle force constant | 20 kJ·mol⁻¹·rad⁻² | |
| Position restraint force | 1000 kJ·mol⁻¹·nm⁻² | Applied under `-DPOSRES` |

**Charge assignment.** Charges are assigned by residue name, with the
protonation state encoded upstream (step 2):

- +1: ARG, LYS, HIP/HSP
- −1: ASP, GLU
- 0: all other residues, including neutral protonation variants (GLH, ASH)

> **Note:** Glycans are not included in the current model. Bead names starting
> with `G` are reserved for future glycan support and are skipped when building
> bonds and angles.

---

## License

[Add your license here, e.g. MIT.]
