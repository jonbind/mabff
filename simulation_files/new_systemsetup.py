

def load_formulation_json(package_dir):
    """Load formulation.json from package directory"""
    json_path = os.path.join(package_dir, 'formulation.json')
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"formulation.json not found in {package_dir}")

    with open(json_path, 'r') as f:
        return json.load(f)

def copy_required_files(package_dir):
    """Copy all required files from requirements folder to package directory"""
    print("="*80)
    print("COPYING REQUIRED FILES")
    print("="*80)

    if not os.path.exists(REQUIREMENTS_FOLDER):
        raise FileNotFoundError(f"Requirements folder not found: {REQUIREMENTS_FOLDER}")

    # List of required files
    required_files = [
        'water.gro',
        'water.itp',
        'ions.itp',
        'mabff220.itp',
        'HST.pdb',
        'CL.pdb',
        'NA.pdb',
        # MDP files
        'minimization.mdp',
        'collect_viscosity',
        'nvt.mdp',
        'npt.mdp',
        'dynamic.mdp',
        'shear.mdp',
        # Analysis script
        'visco_analyse.py',
        'vdwradii.dat'
    ]

    copied_count = 0
    for filename in required_files:
        source = os.path.join(REQUIREMENTS_FOLDER, filename)
        destination = os.path.join(package_dir, filename)

        if os.path.exists(source):
            shutil.copy2(source, destination)
            print(f"✓ Copied: {filename}")
            copied_count += 1
        else:
            print(f"⚠ Warning: {filename} not found in requirements folder")

    print(f"\n✓ Copied {copied_count}/{len(required_files)} required files")
    print("="*80 + "\n")

def mabconccalc(mabconc, dimlen=35.0):
    """Calculate mAb concentration and number of molecules"""
    print(f"Using mAb concentration: {mabconc} mg/mL")
    print(f"Box edge length: {dimlen} nm")

    # Calculate the volume in mL
    volume = (dimlen**3) * 1e-21

    # Calculate the weight in mg
    weight = mabconc * volume

    # Calculate the number of mAb molecules
    mabnr = weight / (148 * 1.660538921e-18)

    # Output the results
    print(f"Volume: {volume:.6e} mL")
    print(f"Weight: {weight:.6e} mg")
    print(f"Number of mAb molecules: {mabnr:.6f}")
    return dimlen, mabnr, mabconc, volume

def bufferconccalc(volume, conc, dissociates=True):
    """Calculate buffer concentration (non-interactive)"""
    volumel = (volume * 1e-3)
    print(f"Box volume: {volume} mL or {volumel} L")
    print(f"Using target buffer concentration: {conc} M")

    moles = (volumel * conc)
    print(f"Moles of buffer in box: {moles}")
    number_buffer = int(round(moles * 6.02214076e23, 0))

    if dissociates:
        print(f"Placing {number_buffer} cations and {number_buffer} anions for the buffer")
    else:
        print(f"Placing total of {number_buffer} buffer molecules")

    return number_buffer

def conccalc(volume, conc, dissociates=True):
    """Calculate salt concentration (non-interactive)"""
    volumel = (volume * 1e-3)
    print(f"Box volume: {volume} mL or {volumel} L")
    print(f"Using target salt concentration: {conc} M")

    moles = (volumel * conc)
    print(f"Moles of salt in box: {moles}")
    numberm = int(round(moles * 6.02214076e23, 0))

    if dissociates:
        print(f"Placing {numberm} cations and {numberm} anions for the salt")
    else:
        print(f"Placing total of {numberm} salt molecules")

    return numberm

def sysprep(dimlen, mabnr, mabconc, number_buffer, name_molecule,
            pdb_file, itp_file, ph):
    """Prepare the system for MD simulation"""
    conc_str = str(int(mabconc))
    gro_file = pdb_file.replace('.pdb', '.gro')
    top_file = f"top_{conc_str}.top"
    dimension = dimlen
    minimization_mdp = "minimization.mdp"
    equilibration_mdp = "equilibration.mdp"
    dynamic_mdp = "dynamic.mdp"
    water_gro = "water.gro"

    # Define ions for buffer and salt based on buffer type
    # For Histidine buffer: HST+ and CL-
    buffer_ion_pos = "HST"
    buffer_ion_neg = "CL"
    salt_ion_pos = "NA"
    salt_ion_neg = "CL"

    # Function to run GROMACS commands
    def run_gmx(command):
        print(f"Running: {command}")
        # Commands run directly - spack environment should be activated before running this script
        os.system(command)

    # 1. Prepare the initial system
    run_gmx(f"rm *#*")
    run_gmx(f"rm *dd*")
    run_gmx(f"gmx_mpi editconf -f {pdb_file} -d 1.0 -o {gro_file}")
    run_gmx(f"gmx_mpi editconf -f {gro_file} -c -o {gro_file} -bt triclinic -box {dimension} {dimension} {dimension}")
    run_gmx(f"echo {1} | gmx_mpi genrestr -f {gro_file} -o posre_{conc_str}.itp")

    # 2. Insert molecules into the box
    run_gmx(f"gmx_mpi insert-molecules -f {gro_file} -ci {gro_file} -o renamed_cg2.gro -nmol {max(0, int(round(mabnr))-1)} -try 5000 -selrpos atom")

    # 3. Writing of the top file
    top_file_content = f"""#include "mabff220.itp"
#include "ions.itp"
#include "{itp_file}"
#ifdef POSRES
#include "posre_{conc_str}.itp"
#endif
#include "water.itp"


[ system ]
{mabconc} mgmL {name_molecule} pH{ph}
[ molecules ]
{name_molecule} {int(round(mabnr))}
"""
    with open(top_file, "w") as f:
        f.write(top_file_content)

    # 3. First minimization
    run_gmx(f"gmx_mpi grompp -p {top_file} -f {minimization_mdp} -c renamed_cg2.gro -o minimization-vac1_{conc_str}.tpr -maxwarn 2 2> charge.log")
    run_gmx(f"gmx_mpi mdrun -deffnm minimization-vac1_{conc_str} -v -ntomp 8")

    # 4. Extract charge and add counter-ions
    total_charge = 0
    with open("charge.log", "r") as f:
        for line in f:
            if "System has non-zero total charge:" in line:
                match = re.search(r"System has non-zero total charge:\s*(-?\d+\.\d+)", line)
                if match:
                    total_charge = int(round(float(match.group(1))))
                    print(f"Total charge of the system: {total_charge}")
                    break
        else:
            print("Warning: Could not find the total charge in the log file.")

    # Logic to add counter-ions
    neutralized_gro = f"neutralized_{conc_str}.gro"
    if total_charge > 0:
        run_gmx(f"gmx_mpi insert-molecules -f minimization-vac1_{conc_str}.gro -ci {salt_ion_neg}.pdb -o {neutralized_gro} -nmol {total_charge} -try 5000 -selrpos atom")
        with open(top_file, "a") as f:
            f.write(f"{salt_ion_neg} {total_charge}\n")
    elif total_charge < 0:
        num_pos_ions = abs(total_charge)
        run_gmx(f"gmx_mpi insert-molecules -f minimization-vac1_{conc_str}.gro -ci {salt_ion_pos}.pdb -o {neutralized_gro} -nmol {num_pos_ions} -try 5000 -selrpos atom")
        with open(top_file, "a") as f:
            f.write(f"{salt_ion_pos} {num_pos_ions}\n")
    else:  # No charge, just copy the file for the next step
        run_gmx(f"cp minimization-vac1_{conc_str}.gro {neutralized_gro}")

    # 5. Add Buffer Ions
    buffer_added_gro1 = f"buffer_added1_{conc_str}.gro"
    buffer_added_gro2 = f"buffer_added2_{conc_str}.gro"
    run_gmx(f"gmx_mpi insert-molecules -f {neutralized_gro} -ci {buffer_ion_pos}.pdb -o {buffer_added_gro1} -nmol {number_buffer} -try 5000 -selrpos atom")
    run_gmx(f"gmx_mpi insert-molecules -f {buffer_added_gro1} -ci {buffer_ion_neg}.pdb -o {buffer_added_gro2} -nmol {number_buffer} -try 5000 -selrpos atom")
    with open(top_file, "a") as f:
        f.write(f"{buffer_ion_pos} {number_buffer}\n{buffer_ion_neg} {number_buffer}\n")

    # This variable will hold the name of the final GRO file before solvation
    final_ion_gro_file = buffer_added_gro2


    # 7. Second minimization (uses the correct final gro file)
    run_gmx(f"gmx_mpi grompp -p {top_file} -f {minimization_mdp} -c {final_ion_gro_file} -o minimization-vac2_{conc_str}.tpr")
    run_gmx(f"gmx_mpi mdrun -deffnm minimization-vac2_{conc_str} -v -ntomp 8 ")

    # 8. Solvation
    run_gmx(f"gmx_mpi solvate -cp minimization-vac2_{conc_str}.gro -cs {water_gro} -radius 0.44 -o solvated_{conc_str}.gro -p {top_file}")
    run_gmx(f"gmx_mpi grompp -p {top_file} -c solvated_{conc_str}.gro -f {minimization_mdp} -o minimization_{conc_str}.tpr")
    print(f"✓ Preparation completed for concentration {mabconc} mg/mL")




# ===== STEP 1: Equilibration (NVT and NPT) =====
echo "Starting equilibration..."

for conc in {" ".join(conc_list)}; do
    echo "=== Equilibrating ${{conc}} mg/mL ==="
    cd conc_${{conc}}

    # Minim
    gmx_mpi grompp -p top_${{conc}}.top -c solvated_${{conc}}.gro -f minimization.mdp -o minimization_${{conc}}.tpr -r solvated_${{conc}}.gro
    gmx_mpi mdrun -deffnm minimization_${{conc}} -v -ntomp $OMP_NUM_THREADS

    # NVT
    gmx_mpi grompp -p top_${{conc}}.top -c minimization_${{conc}}.gro -f nvt.mdp -o nvt_${{conc}}.tpr -maxwarn 3 -r minimization_${{conc}}.gro
    mpiexec -n 1 gmx_mpi mdrun -deffnm nvt_${{conc}} -v -nb gpu -ntomp $OMP_NUM_THREADS -bonded gpu -update gpu -pin on -pme gpu

    # NPT
    gmx_mpi grompp -p top_${{conc}}.top -c nvt_${{conc}}.gro -f npt.mdp -r nvt_${{conc}}.gro -t nvt_${{conc}}.cpt -o npt_${{conc}}.tpr -maxwarn 3
    mpiexec -n 1 gmx_mpi mdrun -deffnm npt_${{conc}} -v -nb gpu -ntomp $OMP_NUM_THREADS -bonded gpu -update gpu -pin on -pme gpu

    cd ..
done

# ===== STEP 2: Determine run number =====
echo "Determining run number..."

# Function to get next available run number for dynamic/shear files
get_next_run_number() {{
    local dir=$1
    local run_num=1

    # Check if any dynamic_*.xtc files exist
    while [ -f "${{dir}}/dynamic_${{run_num}}.xtc" ] || [ -f "${{dir}}/shear_${{run_num}}.xtc" ]; do
        run_num=$((run_num + 1))
    done

    echo $run_num
}}

# Get run number from first concentration folder (they should all be synced)
RUN_NUM=$(get_next_run_number conc_{conc_list[0]})
echo "Using run number: $RUN_NUM"
echo "Output files will be: dynamic_${{RUN_NUM}}.xtc, shear_${{RUN_NUM}}.xtc"

# ===== STEP 3: Prepare production runs =====
echo "Preparing production runs..."

for conc in {" ".join(conc_list)}; do
    cd conc_${{conc}}
    gmx_mpi grompp -p top_${{conc}}.top -c npt_${{conc}}.gro -f dynamic.mdp -r npt_${{conc}}.gro -t npt_${{conc}}.cpt -o dynamic_${{RUN_NUM}}.tpr -maxwarn 3
    cd ..
done

# ===== STEP 4: Run all 4 simulations simultaneously =====
echo "Starting production simulations..."

"""

    # Build multidir list for all concentrations
    multidir_list = " ".join([f"conc_{c}" for c in conc_list])

    slurm_content += f"""
mpiexec -n {len(conc_list)} gmx_mpi mdrun -multidir {multidir_list} -deffnm dynamic_${{RUN_NUM}} -v -nb gpu -ntomp $OMP_NUM_THREADS -bonded gpu -update gpu -pin on -pme gpu

# ===== STEP 5: Prepare shear runs =====
echo "Preparing shear runs..."

for conc in {" ".join(conc_list)}; do
    cd conc_${{conc}}
    gmx_mpi grompp -p top_${{conc}}.top -c dynamic_${{RUN_NUM}}.gro -f shear.mdp -r dynamic_${{RUN_NUM}}.gro -t dynamic_${{RUN_NUM}}.cpt -o shear_${{RUN_NUM}}.tpr -maxwarn 3
    cd ..
done

# ===== STEP 6: Run all {len(conc_list)} shear simulations simultaneously =====
echo "Starting shear simulations..."

mpiexec -n {len(conc_list)} gmx_mpi mdrun -multidir {multidir_list} -deffnm shear_${{RUN_NUM}} -v -nb gpu -ntomp $OMP_NUM_THREADS -bonded gpu -pin on -pme gpu

echo quit | nvidia-cuda-mps-control

# ===== STEP 7: Analysis =====
echo "Starting viscosity analysis for run ${{RUN_NUM}}..."

for conc in {" ".join(conc_list)}; do
    echo "=== Analyzing ${{conc}} mg/mL (run ${{RUN_NUM}}) ==="
    cd conc_${{conc}}
    echo "32" | gmx_mpi energy -f shear_${{RUN_NUM}}.edr -o pressure_${{RUN_NUM}}.xvg
    python ../visco_analyse.py -i shear_from_xvg ./pressure_${{RUN_NUM}}.xvg ./shear.mdp 30 > results_${{RUN_NUM}}.txt
    cd ..
done

echo "✓ All simulations and analysis completed for run ${{RUN_NUM}}: {', '.join(conc_list)} mg/mL"
"""

    with open(script_name, "w") as file:
        file.write(slurm_content)

    print(f"✓ File '{script_name}' created successfully (all concentrations: {', '.join(conc_list)} mg/mL)")
