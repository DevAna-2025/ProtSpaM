#!/bin/bash
#SBATCH --job-name=protspam_v55
#SBATCH --output=logs_variants_55/slurm_%x_%j.out
#SBATCH --error=logs_variants_55/slurm_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=48:00:00
#SBATCH --mem=64G

module load cesga/2025
module load gcc
module load openmpi/5.0.9

mkdir -p logs_variants_55 results_variants_55

BIN_NAME="$1"
PHASE3_DIST="$2"
PHASE4_MODE="$3"
CALC_VARIANT="$4"
BATCH_SIZE="${5:-NA}"

if [ -z "$BIN_NAME" ] || [ -z "$PHASE3_DIST" ] || [ -z "$PHASE4_MODE" ] || [ -z "$CALC_VARIANT" ]; then
    echo "Usage:"
    echo "  sbatch run_variant_55.sbatch <bin_name> <phase3_dist> <phase4_mode> <calc_variant> [batch_size]"
    echo
    echo "Example:"
    echo "  sbatch run_variant_55.sbatch protspam_balanced_pipeline_calcbase balanced pipeline calcbase"
    echo "  sbatch run_variant_55.sbatch protspam_block_pipeline_metacache_batch_calcbase block pipeline_metacache_batch2 calcbase 2"
    exit 2
fi

if [ "$BATCH_SIZE" != "NA" ]; then
    export PROTSPAM_BATCH_SIZE="$BATCH_SIZE"
fi

BIN="./bin/Debug/${BIN_NAME}"
FILELIST="filelist_55"
PATTERNS="patterns_clean.txt"
SPECIES="55"
REPS=3
NPS="1 2 4 8 16 32 64"

VARIANT="${PHASE3_DIST}_${PHASE4_MODE}_${CALC_VARIANT}"
SUMMARY="logs_variants_55/resumen_${VARIANT}_${SLURM_JOB_ID}.tsv"

if ! command -v mpirun >/dev/null 2>&1; then
    echo "ERROR: mpirun not found. Check MPI module loading."
    exit 127
fi

if [ ! -x "$BIN" ]; then
    echo "ERROR: executable not found or not executable: $BIN"
    echo "Run first:"
    echo "  make clean"
    echo "  make variants"
    exit 127
fi

if [ ! -f "$FILELIST" ]; then
    echo "ERROR: filelist not found: $FILELIST"
    exit 2
fi

if [ ! -f "$PATTERNS" ]; then
    echo "ERROR: patterns file not found: $PATTERNS"
    exit 2
fi

echo -e "variant\tbin_name\tphase3_dist\tphase4_mode\tcalc_variant\tbatch_size\tjob_id\tnodes\tntasks_alloc\tspecies\tnp\trep\tstatus\texit_code\tt_phase3_s\tt_phase4_s\tt_total_program_s\telapsed_seconds\tmem_per_node\tnodelist\tlog_file\toutput_file" > "$SUMMARY"

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
    np="$1"
    rep="$2"

    name="${VARIANT}_${SPECIES}sp_np${np}_rep${rep}"
    log="logs_variants_55/${name}_${SLURM_JOB_ID}.log"
    out="results_variants_55/DMat_${name}_${SLURM_JOB_ID}"

    echo "==================================================" | tee "$log"
    echo "Running $name" | tee -a "$log"
    echo "Binary: $BIN" | tee -a "$log"
    echo "Variant: $VARIANT" | tee -a "$log"
    echo "Phase 3 distribution: $PHASE3_DIST" | tee -a "$log"
    echo "Phase 4 mode: $PHASE4_MODE" | tee -a "$log"
    echo "calc_matches: $CALC_VARIANT" | tee -a "$log"
    echo "Batch size: $BATCH_SIZE" | tee -a "$log"
    echo "Filelist: $FILELIST" | tee -a "$log"
    echo "Species: $SPECIES" | tee -a "$log"
    echo "NP: $np" | tee -a "$log"
    echo "Rep: $rep" | tee -a "$log"
    echo "SLURM_JOB_ID: $SLURM_JOB_ID" | tee -a "$log"
    echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES" | tee -a "$log"
    echo "SLURM_NTASKS: $SLURM_NTASKS" | tee -a "$log"
    echo "SLURM_MEM_PER_NODE: $SLURM_MEM_PER_NODE" | tee -a "$log"
    echo "SLURM_NODELIST: $SLURM_NODELIST" | tee -a "$log"
    echo "Start: $(date)" | tee -a "$log"
    echo "==================================================" | tee -a "$log"

    start_sec=$(date +%s)

    mpirun -np "$np" "$BIN" \
        -l "$FILELIST" \
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

    echo -e "${VARIANT}\t${BIN_NAME}\t${PHASE3_DIST}\t${PHASE4_MODE}\t${CALC_VARIANT}\t${BATCH_SIZE}\t${SLURM_JOB_ID}\t${SLURM_JOB_NUM_NODES}\t${SLURM_NTASKS}\t${SPECIES}\t${np}\t${rep}\t${status}\t${code}\t${t_phase3}\t${t_phase4}\t${t_total}\t${elapsed}\t${SLURM_MEM_PER_NODE}\t${SLURM_NODELIST}\t${log}\t${out}" >> "$SUMMARY"
}

for np in $NPS; do
    for rep in $(seq 1 "$REPS"); do
        run_test "$np" "$rep"
    done
done

echo "=================================================="
echo "Resumen guardado en: $SUMMARY"
echo "=================================================="
cat "$SUMMARY"
