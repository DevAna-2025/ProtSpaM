# ProtSpaM-MPI

MPI extension of **Prot-SpaM** developed as part of a Master's Thesis in High Performance Computing.

This branch contains the MPI implementation used for the Phase 4 experiments. It keeps the Phase 3 Option A approach for input handling and adds a memory-aware parallel Phase 4 based on pattern-level streaming and optimized match computation.

## Original Project

This work is based on:

> Leimeister, C. A., Schellhorn, J., Schoebel, S., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).
> *Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.*
> bioRxiv, 306142.

Original repository:

https://github.com/jschellh/ProtSpaM

## Current MPI Implementation

The implementation parallelizes the most expensive stages of Prot-SpaM:

- **Phase 3: spaced-word generation**
  - Species are distributed across MPI ranks.
  - The distribution is balanced using sequence length, instead of assigning contiguous blocks of species.
  - Each rank computes spaced-words for its local species.

- **Phase 4: match computation**
  - Distance matrix pairs are computed in parallel.
  - Spaced-words are processed in a pattern-level pipeline to reduce peak memory usage.
  - Remote spaced-words are streamed pattern by pattern instead of replicating all species data on every rank.
  - Equal-key run lengths are precomputed to avoid repeated `multiMatch` scans inside the match kernel.
  - Don't-care positions are precomputed per pattern and reused during scoring.

This version is referred to in the benchmark scripts as:

```text
pipeline_runlen
```

## Requirements

- C++ compiler with C++11 support
- MPI compiler wrapper, such as `mpicxx`
- GNU Make

On FinisTerrae III, the experiments were run with:

```bash
module load cesga/2025
module load gcc
module load openmpi/5.0.9
```

## Compilation

```bash
make clean
make
```

The executable is generated at:

```bash
./bin/Debug/protspam
```

## Preparing the Data

Create a `data/` directory and place the FASTA/proteome files according to the paths used in the filelists.

The biological datasets are not included in this repository. The files listed in:

- `filelist_10`
- `filelist_20`
- `filelist_30`

must exist before running the experiments.

Example:

```bash
mkdir -p data
cp /path/to/proteomes/*.faa data/
```

If a filelist contains:

```text
data/species1.faa
data/species2.faa
```

then those files must exist at those paths.

## Running Manually

```bash
mpirun -np <processes> ./bin/Debug/protspam \
    -l <filelist> \
    -p patterns_clean.txt \
    -o <output_matrix>
```

Example:

```bash
mpirun -np 4 ./bin/Debug/protspam \
    -l filelist_20 \
    -p patterns_clean.txt \
    -o DMat_20sp_np4
```

The program reports separate timings for:

```text
Tiempo spaced-words
Tiempo matches
Tiempo total
```

In MPI executions, Phase 3 and Phase 4 are internally interleaved by pattern to reduce memory usage, but the accumulated times are still reported separately.

## Benchmark Scripts

Two SLURM scripts are included for Phase 4 benchmarking.

### Single-node benchmark

```bash
sbatch run_phase4_single.sbatch
```

This script evaluates:

- Species sets: 10, 20 and 30
- MPI processes: 1, 2, 4, 8, 16 and 32
- Repetitions: 5
- Nodes: 1

Outputs:

```text
logs_phase4_pipeline/
results_phase4_pipeline/
logs_phase4_pipeline/resumen_single_<job_id>.tsv
```

### Multinode benchmark

```bash
sbatch run_phase4_multinode.sbatch
```

This script evaluates:

- Species sets: 20 and 30
- MPI processes: 8 and 16
- Repetitions: 3
- Nodes: 2

Outputs:

```text
logs_phase4_multinode/
results_phase4_multinode/
logs_phase4_multinode/resumen_multi_<job_id>.tsv
```

The multinode benchmark is intended as an exploratory scalability experiment. The main performance analysis should use the single-node benchmark, while multinode results are useful for studying memory pressure and communication overhead.

## Output Files

The main output is a distance matrix file:

```text
DMat_*
```

Benchmark scripts also generate:

- one execution log per run;
- one output matrix per run;
- one TSV summary containing status, exit code, Phase 3 time, Phase 4 time, total program time and elapsed wall-clock time.

Generated logs, matrices and binaries are ignored by Git.

## Repository Structure

```text
.
|-- include/
|-- src/
|-- main.cpp
|-- Makefile
|-- filelist_10
|-- filelist_20
|-- filelist_30
|-- patterns_clean.txt
|-- run_phase4_single.sbatch
|-- run_phase4_multinode.sbatch
|-- run_opcionA.sbatch
|-- README.md
|-- COPYING
```

`run_opcionA.sbatch` is kept as a historical Phase 3 Option A benchmark script. The Phase 4 results reported for this version should use `run_phase4_single.sbatch` and `run_phase4_multinode.sbatch`.

## Notes on the Final Phase 4 Version

The final Phase 4 implementation was designed after observing memory failures in earlier prototypes. The previous approaches either gathered too much data on rank 0 or replicated too much spaced-word data across MPI ranks. The current pipeline version reduces peak memory by consuming one pattern at a time and discarding temporary spaced-word data before moving to the next pattern.

This improves robustness for larger species sets while preserving the output distance matrix, as verified by comparing MPI outputs against one-process reference outputs.
