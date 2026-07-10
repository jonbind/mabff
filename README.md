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
