#!/bin/bash
#SBATCH --job-name=protspam_multinodo
#SBATCH --output=logs_multinodo/slurm_%j.out
#SBATCH --error=logs_multinodo/slurm_%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH --mem=64G

module load cesga/2025
module load gcc
module load openmpi/5.0.9

mkdir -p logs_multinodo results_multinodo

# BENCHMARK MULTINODO (Opción A: escalado por número de nodos).
# Densidad fija de 32 procesos por nodo; se varía el número de nodos:
#   1 nodo  ->  32 procesos
#   2 nodos ->  64 procesos
#   4 nodos -> 128 procesos
#   8 nodos -> 256 procesos
# Se reserva el máximo (8 nodos) y dentro del job se lanza cada configuración
# pidiendo el número de nodos/procesos correspondiente con mpirun.
# Dataset: filelist_64. Las tres variantes. 3 repeticiones.
FILELIST="filelist_64"
PATTERNS="patterns_clean.txt"
PPN=32                      # procesos por nodo
NODES_LIST="1 2 4 8"        # número de nodos a evaluar
REPS=3

declare -A BINS
BINS["metacache"]="./bin/Debug/protspam_block_pipeline_metacache_calcbase"
BINS["isend"]="./bin/Debug/protspam_block_pipeline_metacache_isend_calcbase"
BINS["isend_opt"]="./bin/Debug/protspam_block_pipeline_metacache_isend_calcopt"
VARIANTS="metacache isend isend_opt"

SUMMARY="logs_multinodo/resumen_multinodo_${SLURM_JOB_ID}.tsv"

if ! command -v mpirun >/dev/null 2>&1; then echo "ERROR: mpirun not found."; exit 127; fi
for v in $VARIANTS; do
    if [ ! -x "${BINS[$v]}" ]; then echo "ERROR: no existe ${BINS[$v]}"; exit 127; fi
done

echo -e "variante\tjob_id\tnodos\tppn\tnp_total\trep\tstatus\texit_code\tt_phase3_s\tt_phase4_s\tt_total_program_s\telapsed_seconds\tlog_file\toutput_file" > "$SUMMARY"

extract_time() { value=$(grep "$1" "$2" | tail -n 1 | awk '{print $(NF-1)}'); if [ -z "$value" ]; then echo "NA"; else echo "$value"; fi; }

run_test() {
    variant="$1"; nodes="$2"; rep="$3"; bin="${BINS[$variant]}"
    np=$(( nodes * PPN ))
    name="${variant}_64sp_${nodes}nodos_np${np}_rep${rep}"
    log="logs_multinodo/${name}_${SLURM_JOB_ID}.log"
    out="results_multinodo/DMat_${name}_${SLURM_JOB_ID}"
    echo "=== $name ===" | tee "$log"; echo "Start: $(date)" | tee -a "$log"
    start=$(date +%s)
    # -N fija el número de nodos, -np el total de procesos; --map-by node reparte
    mpirun -np "$np" --map-by ppr:${PPN}:node "$bin" -l "$FILELIST" -p "$PATTERNS" -o "$out" >> "$log" 2>&1
    code=$?; end=$(date +%s); elapsed=$((end-start))
    echo "End: $(date)" >> "$log"; echo "Exit code: $code" >> "$log"; echo "Elapsed seconds: $elapsed" >> "$log"
    [ "$code" -eq 0 ] && status="OK" || status="FAIL"
    t3=$(extract_time "Tiempo spaced-words:" "$log"); t4=$(extract_time "Tiempo matches:" "$log"); tt=$(extract_time "Tiempo total:" "$log")
    echo -e "${variant}\t${SLURM_JOB_ID}\t${nodes}\t${PPN}\t${np}\t${rep}\t${status}\t${code}\t${t3}\t${t4}\t${tt}\t${elapsed}\t${log}\t${out}" >> "$SUMMARY"
}

for variant in $VARIANTS; do
    for nodes in $NODES_LIST; do
        for rep in $(seq 1 "$REPS"); do run_test "$variant" "$nodes" "$rep"; done
    done
done

echo "=== Resumen: $SUMMARY ==="; cat "$SUMMARY"
