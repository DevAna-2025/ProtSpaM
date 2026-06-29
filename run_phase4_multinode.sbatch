#!/bin/bash
#SBATCH --job-name=protspam_phase4_multi
#SBATCH --output=logs_phase4_multinode/slurm_multi_%j.out
#SBATCH --error=logs_phase4_multinode/slurm_multi_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G

module load cesga/2025
module load gcc
module load openmpi/5.0.9

mkdir -p logs_phase4_multinode results_phase4_multinode

BIN="./bin/Debug/protspam"
PATTERNS="patterns_clean.txt"

VERSION="pipeline_runlen"
REPS=3
NPS="8 16"

SUMMARY="logs_phase4_multinode/resumen_multi_${SLURM_JOB_ID}.tsv"

if ! command -v mpirun >/dev/null 2>&1; then
    echo "ERROR: mpirun not found. Check MPI module loading."
    exit 127
fi

if [ ! -x "$BIN" ]; then
    echo "ERROR: executable not found or not executable: $BIN"
    echo "Run: make clean && make"
    exit 127
fi

echo -e "version\tmode\tjob_id\tnodes\tntasks_alloc\tntasks_per_node\tspecies\tnp\trep\tstatus\texit_code\tt_phase3_s\tt_phase4_s\tt_total_program_s\telapsed_seconds\tmem_per_node\tnodelist\tlog_file\toutput_file" > "$SUMMARY"

extract_time() {
    pattern="$1"
    log_file="$2"

    value=$(grep "$pattern" "$log_file" | tail -n 1 | awk '{print $(NF-1)}')

    if [ -z "$value" ]; then
        echo "NA"
    else
        echo "$value"
    fi
}

run_test() {
    species="$1"
    np="$2"
    rep="$3"
    filelist="$4"

    name="${VERSION}_multi_${species}sp_np${np}_rep${rep}"
    log="logs_phase4_multinode/${name}_${SLURM_JOB_ID}.log"
    out="results_phase4_multinode/DMat_${name}_${SLURM_JOB_ID}"

    echo "==================================================" | tee "$log"
    echo "Running $name" | tee -a "$log"
    echo "Version: $VERSION" | tee -a "$log"
    echo "Mode: multinode" | tee -a "$log"
    echo "Filelist: $filelist" | tee -a "$log"
    echo "Species label: $species" | tee -a "$log"
    echo "NP: $np" | tee -a "$log"
    echo "Rep: $rep" | tee -a "$log"
    echo "SLURM_JOB_ID: $SLURM_JOB_ID" | tee -a "$log"
    echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES" | tee -a "$log"
    echo "SLURM_NTASKS: $SLURM_NTASKS" | tee -a "$log"
    echo "SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE" | tee -a "$log"
    echo "SLURM_MEM_PER_NODE: $SLURM_MEM_PER_NODE" | tee -a "$log"
    echo "SLURM_NODELIST: $SLURM_NODELIST" | tee -a "$log"
    echo "Start: $(date)" | tee -a "$log"
    echo "==================================================" | tee -a "$log"

    start_sec=$(date +%s)

    mpirun -np "$np" "$BIN" \
        -l "$filelist" \
        -p "$PATTERNS" \
        -o "$out" \
        >> "$log" 2>&1

    code=$?

    end_sec=$(date +%s)
    elapsed=$((end_sec - start_sec))

    echo "End: $(date)" >> "$log"
    echo "Elapsed seconds: $elapsed" >> "$log"
    echo "Exit code: $code" >> "$log"

    if [ "$code" -eq 0 ]; then
        status="OK"
    else
        status="FAIL"
    fi

    t_phase3=$(extract_time "Tiempo spaced-words:" "$log")
    t_phase4=$(extract_time "Tiempo matches:" "$log")
    t_total=$(extract_time "Tiempo total:" "$log")

    echo -e "${VERSION}\tmultinode\t${SLURM_JOB_ID}\t${SLURM_JOB_NUM_NODES}\t${SLURM_NTASKS}\t${SLURM_NTASKS_PER_NODE}\t${species}\t${np}\t${rep}\t${status}\t${code}\t${t_phase3}\t${t_phase4}\t${t_total}\t${elapsed}\t${SLURM_MEM_PER_NODE}\t${SLURM_NODELIST}\t${log}\t${out}" >> "$SUMMARY"
}

run_species_set() {
    species="$1"
    filelist="$2"

    for np in $NPS; do
        for rep in $(seq 1 "$REPS"); do
            run_test "$species" "$np" "$rep" "$filelist"
        done
    done
}

run_species_set 20 filelist_20
run_species_set 30 filelist_30

echo "=================================================="
echo "Resumen guardado en: $SUMMARY"
echo "=================================================="
cat "$SUMMARY"
