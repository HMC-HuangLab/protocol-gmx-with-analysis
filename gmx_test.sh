#!/bin/bash
#run test.sh 500 pdb_dir LIG 1 /mdp_dir/


SIMULATIONTIME=$1
PDB_DIR=$2
LIGNAME=$3
GPU_ID=$4
MDP_DIR=$5 # 新增的参数，用于指定 MDP 文件的目录

#---------  SIMU SETUP  -----------
echo "SIMULATIONTIME: $SIMULATIONTIME"
echo "PDB_DIR: $PDB_DIR"
echo "LIGNAME: $LIGNAME"
echo "GPU_ID: $GPU_ID"
echo "MDP_DIR: $MDP_DIR" # 打印出 MDP 目录信息


#---------  SIMU SETUP  -----------
BOXSIZE=1.5 # cubic simulation box size
BOXTYPE=cubic # Box type
NT=16 # Number of cores
WATER=tip3p # Water type
NUMBEROFREPLICAS=1 # Number of replicas
FF=amber99sb-ildn # Force field

#---------  HPC SETUP  -----------
MPI="" # If you have to submit jobs with MPI software like "mpirun -np 10", add the command here
GMX=gmx # GMX command (can be "$GMX_mpi" sometimes, just change it here)
# GPU configuration
GPU0="-gpu_id $GPU_ID -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -npme 1 -pme gpu -pmefft gpu -pin on -pinstride 0 -nstlist 100 -pinoffset 0"
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_GPU_DD_COMMS=true
MDRUN_CPU="$GMX mdrun -nt $NT"
MDRUN_GPU="$GMX mdrun $GPU0"

# Log file for recording progress and errors
LOGFILE="simulation.log"

log() {
    echo "$(date): $*" >> $LOGFILE
}

# Function to handle errors and cleanup
error_handler() {
    log "An error occurred in $FILE. Skipping to the next PDB file."
    cd ../..
    return 1
}

# Function to process each PDB file
process_pdb_file() {
    local PDB_FILE=$1
    local FILE=$(basename "$PDB_FILE" .pdb)

    mkdir -p "$FILE"
    cp "$PDB_FILE" "$FILE"/
    if [ ! -d "$MDP_DIR" ]; then
        log "MDP directory not found: $MDP_DIR!"
        cd ..
        return 1
    fi
    cp -r "$MDP_DIR" "$FILE"/  
    cd "$FILE"

    # Setup simulation time
    if [ ! -z "$SIMULATIONTIME" ]; then
        python <<EOF
import re
restep = re.compile("nsteps *= *(\d*)")
redt = re.compile("dt *= *(\d*.\d*)")
dt = 0
simulationtime = float($SIMULATIONTIME) *1000 # Time in ps
outputLines = []
with open("mdp/md_prod.mdp", 'r') as f:
    mdp = f.readlines()
    # Find the timestep first
    for line in mdp: 
        dtmatch = redt.match(line)
        if dtmatch:
            dt = float(dtmatch.group(1))
            break
    for line in mdp:
        stepmatch = restep.match(line)
        if stepmatch and float(dt) > 0:
            nsteps = int(simulationtime) / dt
            line = "nsteps            = {}        ; {} * {} = {} ps or {} ns\n".format(int(nsteps), dt, nsteps, dt * nsteps, simulationtime / 1000)
        outputLines.append(line)
    with open("mdp/md_prod.mdp", 'w') as f:
        for line in outputLines:
            f.write(line)
EOF
    fi

    if [ ! -z "$LIGNAME" ]; then
        mkdir -p param_"$FILE"
        cp "$PDB_FILE" param_"$FILE"/
        cd param_"$FILE"

        grep 'ATOM  ' "$FILE.pdb" --color=none > receptor.pdb

        # Extract ligand and connect
        python <<EOF
ligand_atom = []
keepLine = []
with open("$FILE.pdb", "r") as file:
    lines = file.readlines()
    for line in lines:
        if '$LIGNAME' in line[17:20]:
            line = line[:17] + "LIG" + line[20:]
            keepLine.append(line)
            ligand_atom.append(int(line[6:11]))
        elif "CONECT" in line[0:6]:
            idx = [int(x) for x in line.split()[1:]]
            if any(id in idx for id in ligand_atom):
                keepLine.append(line)
with open("ligand.pdb", "w") as file:
    for line in keepLine:
        file.write(line)
EOF

        if [ ! -d "ligand" ]; then
            mkdir ligand
        fi
        if [ ! -d "ligand/ligand.acpype" ]; then
            obabel -ipdb ligand.pdb -omol2 -h > ligand.mol2
            acpype -i ligand.mol2
            mv ligand* ligand/
        fi

        if [ ! -d "receptor" ]; then
            mkdir receptor
        fi
        mv receptor.pdb receptor/
        cd receptor

        $GMX pdb2gmx -f receptor.pdb -o receptor_GMX.pdb -water $WATER -ignh -ff $FF
        if [ $? -ne 0 ]; then
            log "pdb2gmx failed for $FILE"
            cd ../..
            return 1
        fi

        cd ../../
        cp param_"$FILE"/receptor/*.itp param_"$FILE"/receptor/topol.top .
        cp param_"$FILE"/ligand/ligand.acpype/ligand_GMX.itp ligand.itp
        grep -h ATOM param_"$FILE"/receptor/receptor_GMX.pdb param_"$FILE"/ligand/ligand.acpype/ligand_NEW.pdb > complex.pdb

        cp topol.top topol.bak
        sed '/forcefield\.itp\"/a\
#include "ligand.itp"
' topol.top > topol2.top
        mv topol2.top topol.top
        echo "ligand   1" >> topol.top

        ndx=$($GMX make_ndx -f param_"$FILE"/ligand/ligand.acpype/ligand_NEW.pdb -o lig_noh.ndx <<EOF
r LIG & !a H*
name 3 LIG-H
q
EOF
)
        cp param_"$FILE"/ligand/ligand.acpype/posre_ligand.itp .

        echo '
; Include Position restraint file
#ifdef POSRES
#include "posre_ligand.itp"
#endif' >> ligand.itp

        $GMX editconf -f complex.pdb -o complex_newbox.gro -d $BOXSIZE -bt $BOXTYPE
        if [ $? -ne 0 ]; then
            log "editconf failed for $FILE"
            cd ../..
            return 1
        fi
        PDB=complex
    else
        PDB=$FILE
        $GMX pdb2gmx -f "$PDB.pdb" -o "$PDB""_processed.gro" -water $WATER -ignh -ff $FF
        if [ $? -ne 0 ]; then
            log "pdb2gmx failed for $FILE"
            cd ..
            return 1
        fi
        $GMX editconf -f "$PDB""_processed.gro" -o "$PDB""_newbox.gro" -d $BOXSIZE -bt $BOXTYPE
        if [ $? -ne 0 ]; then
            log "editconf failed for $FILE"
            cd ..
            return 1
        fi
    fi

    $GMX solvate -cp "$PDB""_newbox.gro" -cs spc216.gro -o "$PDB""_solv.gro" -p topol.top
    if [ $? -ne 0 ]; then
        log "solvate failed for $FILE"
        cd ..
        return 1
    fi

    $GMX grompp -f mdp/ions.mdp -c "$PDB""_solv.gro" -p topol.top -o ions.tpr --maxwarn 1
    if [ $? -ne 0 ]; then
        log "grompp failed for ions for $FILE"
        cd ..
        return 1
    fi
    echo "SOL" | $GMX genion -s ions.tpr -o "$PDB""_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral
    if [ $? -ne 0 ]; then
        log "genion failed for $FILE"
        cd ..
        return 1
    fi

    if [ ! -z "$LIGNAME" ]; then
        ndx=$($GMX make_ndx -f complex_solv_ions.gro -o index.ndx <<EOF
1 | r LIG
r SOL | r CL | r NA
q
EOF
)
        if [ $? -ne 0 ]; then
            log "make_ndx failed for $FILE"
            cd ..
            return 1
        fi

        python <<EOF
import re
with open('index.ndx', 'r') as file:
    content = file.read()
matches = re.findall(r'\[ \w+ \]', content)
if matches:
    content = content.replace(matches[-1], '[ Water_Ions ]')
    content = content.replace(matches[-2], '[ Protein_Ligand ]')
    with open('index.ndx', 'w') as file:
        file.write(content)
EOF

        sed -i 's/Protein Non-Protein/Protein_Ligand Water_Ions/g' mdp/nvt_300.mdp
        sed -i 's/Protein Non-Protein/Protein_Ligand Water_Ions/g' mdp/npt.mdp
        sed -i 's/Protein Non-Protein/Protein_Ligand Water_Ions/g' mdp/md_prod.mdp

        INDEX="-n index.ndx"
    else
        INDEX=""
    fi

    for ((i = 0; i < $NUMBEROFREPLICAS; i++)); do
        echo ">>>>> replica_$((i + 1))_$FILE"
        mkdir -p "replica_$((i + 1))_$FILE/graph"
        mkdir -p "replica_$((i + 1))_$FILE/gro"
        mkdir -p "replica_$((i + 1))_$FILE/mdp"
        mkdir -p "replica_$((i + 1))_$FILE/results/mini"
        mkdir -p "replica_$((i + 1))_$FILE/results/nvt"
        mkdir -p "replica_$((i + 1))_$FILE/results/npt"
        mkdir -p "replica_$((i + 1))_$FILE/results/prod"

        cd "replica_$((i + 1))_$FILE"
        cp -R ../mdp .
        cp ../"$PDB""_solv_ions.gro" .
        cp ../topol.top .
        cp ../*.itp .
        cp ../index.ndx . 2> /dev/null

        log "Running energy minimization for replica_$((i + 1))_$FILE"
        $GMX grompp -f mdp/em.mdp -c "$PDB""_solv_ions.gro" -p topol.top -o em.tpr $INDEX
        if ! $MPI $MDRUN_CPU -v -deffnm em; then
            log "Energy minimization failed for replica_$((i + 1))_$FILE"
            cd ../..
            return 1
        fi

        mv em* results/mini/
        mv mdout.mdp results/mini/
        mv *.gro gro/
        echo "11 0" | $GMX energy -f results/mini/em.edr -o graph/mini_"$PDB"_pot.xvg

        # Check energy minimization results
        energy_check=$(gmx energy -f results/mini/em.edr <<< "10 0" | grep Potential | awk '{print $2}')
        if (( $(echo "$energy_check > 1.0e+06" |bc -l) )); then
            log "Energy minimization did not converge. Potential energy is too high: $energy_check"
            cd ../..
            return 1
        fi

        log "Running NVT equilibration for replica_$((i + 1))_$FILE"
        $GMX grompp -f mdp/nvt_300.mdp -c results/mini/em.gro -r results/mini/em.gro -p topol.top -o nvt_300.tpr -maxwarn 2 $INDEX
        if ! $MPI $MDRUN_GPU -deffnm nvt_300 -v; then
            log "NVT equilibration failed for replica_$((i + 1))_$FILE"
            cd ../..
            return 1
        fi

        mv nvt* results/nvt/ 2> /dev/null
        mv mdout.mdp results/nvt/
        echo "16 0" | $GMX energy -f results/nvt/nvt_300.edr -o graph/temperature_nvt_300.xvg

        log "Running NPT equilibration for replica_$((i + 1))_$FILE"
        $GMX grompp -f mdp/npt.mdp -c results/nvt/nvt_300.gro -r results/nvt/nvt_300.gro -t results/nvt/nvt_300.cpt -p topol.top -o npt_ab.tpr -maxwarn 2 $INDEX
        if ! $MPI $MDRUN_GPU -deffnm npt_ab -v; then
            log "NPT equilibration failed for replica_$((i + 1))_$FILE"
            cd ../..
            return 1
        fi

        mv npt* results/npt/ 2> /dev/null
        mv mdout.mdp results/npt_ab/
        echo "17 0" | $GMX energy -f results/npt/npt_ab.edr -o graph/npt_"$PDB"_pressure.xvg
        echo "22 0" | $GMX energy -f results/npt/npt_ab.edr -o graph/npt_"$PDB"_volume.xvg

        log "Running production MD for replica_$((i + 1))_$FILE"
        $GMX grompp -f mdp/md_prod.mdp -c results/npt/npt_ab.gro -t results/npt/npt_ab.cpt -p topol.top -o "md_"$PDB"_prod.tpr" -maxwarn 2 $INDEX
        if ! $MPI $MDRUN_GPU -deffnm "md_"$PDB"_prod" -v; then
            log "Production MD failed for replica_$((i + 1))_$FILE"
            cd ../..
            return 1
        fi

        mv md_* results/prod 2> /dev/null
        mv mdout.mdp results/prod/
        echo "backbone backbone" | $GMX rms -s "results/prod/md_"$PDB"_prod.tpr" -f "results/prod/md_"$PDB"_prod.trr" -o graph/prod_"$PDB"_rmsd.xvg -tu ns

        cd results/prod
        echo "Protein System" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_prod.trr" -o "md_"$PDB"_clean_temp.xtc" -pbc nojump -ur compact -center
        echo "Protein System" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_full.xtc" -fit rot+trans
        echo "Protein non-Water" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_nowat.xtc" -fit rot+trans
        echo "Protein non-Water" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_nowat.pdb" -pbc nojump -ur compact -center -b 0 -e 0

        echo "Protein Protein" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_protein_LF.pdb" -pbc nojump -ur compact -center -dump 9999999999999999
        rm "md_"$PDB"_clean_temp.pdb"
        echo "non-Water" | $GMX convert-tpr -s "md_"$PDB"_prod.tpr" -o tpr_nowat.tpr
        echo "Protein" | $GMX -s tpr_nowat.tpr -f "md_"$PDB"_clean_nowat.xtc" -ol "md_"$PDB"_clean_nowat_filtered.xtc" -all -fit

        cd ../../
        echo "backbone backbone" | $GMX rms -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_rmsd_ca.xvg" -tu ns
        echo "protein protein" | $GMX rms -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_rmsd_all.xvg" -tu ns
        echo "backbone" | $GMX gyrate -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_gyrate.xvg"
        echo "backbone" | $GMX rmsf -s "results/prod/tpr_nowat.tpr" -f "results/prod/md_"$PDB"_clean_nowat.xtc" -o "graph/prod_"$PDB"_rmsf_ref.xvg" -res

        cd ../
    done

    cd ../
}

# Main loop to process all PDB files
for PDB_FILE in $PDB_DIR/*.pdb; do
    process_pdb_file "$PDB_FILE"
    if [ $? -ne 0 ]; then
        log "Skipping $PDB_FILE due to error"
        continue
    fi
done
