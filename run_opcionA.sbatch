#!/bin/bash
#SBATCH --job-name=protspam_A_benchmark
#SBATCH --output=logs/slurm_A_benchmark_%j.out
#SBATCH --error=logs/slurm_A_benchmark_%j.err
#SBATCH --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem=64G

mkdir -p logs results

BIN="./bin/Debug/protspam"
PATTERNS="patterns_clean.txt"

REPS=5
NPS="1 2 4 8 16 32"

SUMMARY="logs/resumen_opcionA_${SLURM_JOB_ID}.tsv"

echo -e "species\tnp\trep\tstatus\texit_code\telapsed_seconds\tlog_file\toutput_file" > "$SUMMARY"

run_test() {
    species="$1"
    np="$2"
    rep="$3"
    filelist="$4"

    name="A_${species}sp_np${np}_rep${rep}"
    log="logs/${name}_${SLURM_JOB_ID}.log"
    out="results/DMat_${name}_${SLURM_JOB_ID}"

    echo "==================================================" | tee "$log"
    echo "Running $name" | tee -a "$log"
    echo "Filelist: $filelist" | tee -a "$log"
    echo "Species label: $species" | tee -a "$log"
    echo "NP: $np" | tee -a "$log"
    echo "Rep: $rep" | tee -a "$log"
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

    echo -e "${species}\t${np}\t${rep}\t${status}\t${code}\t${elapsed}\t${log}\t${out}" >> "$SUMMARY"
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

run_species_set 10 filelist_10
run_species_set 20 filelist_20
run_species_set 30 filelist_30

echo "=================================================="
echo "Resumen guardado en: $SUMMARY"
echo "=================================================="
cat "$SUMMARY"
