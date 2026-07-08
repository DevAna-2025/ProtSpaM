# ProtSpaM-MPI

MPI extension of **Prot-SpaM** developed as part of a Master's Thesis in High Performance Computing.

This implementation parallelizes the most expensive stages of Prot-SpaM while reducing memory consumption during match computation through a pattern-level streaming pipeline.

## Original Project

Based on:

> Leimeister, C. A., Schellhorn, J., Schoebel, S., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).  
> *Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.*

Original repository:

https://github.com/jschellh/ProtSpaM

---

## MPI Implementation

### Phase 3 – Spaced-word generation

- Species are distributed across MPI ranks.
- Load balancing uses sequence length.
- Each rank computes spaced-words for its assigned species.

### Phase 4 – Match computation

- Species pairs are processed in parallel.
- Spaced-words are streamed one pattern at a time to reduce memory usage.
- Remote spaced-words are exchanged on demand instead of replicating the complete dataset.
- Equal-key run lengths and don't-care positions are precomputed to accelerate matching.

This implementation corresponds to the benchmark version:

```text
pipeline_runlen
```

---

## Requirements

- C++11 compiler
- MPI implementation (e.g. OpenMPI)
- GNU Make

Example environment on FinisTerrae III:

```bash
module load cesga/2025
module load gcc
module load openmpi/5.0.9
```

---

## Compilation

```bash
make clean
make
```

Executable:

```text
./bin/Debug/protspam
```

---

## Preparing the Data

Create the data directory and copy the proteome files referenced by the filelists.

```bash
mkdir -p data
cp /path/to/proteomes/*.faa data/
```

Datasets are **not included** in this repository.

The experiments used the following filelists:

- `filelist_10`
- `filelist_20`
- `filelist_30`
- `filelist_50`
- `filelist_55`

---

## Running

Example:

```bash
mpirun -np 4 ./bin/Debug/protspam \
    -l filelist_20 \
    -p patterns_clean.txt \
    -o DMat_20sp_np4
```

The program reports:

```text
Tiempo spaced-words
Tiempo matches
Tiempo total
```

---

## Benchmarks

### Single-node

Create the output directories:

```bash
mkdir -p logs_phase4_single_64g results_phase4_single_64g
```

Run:

```bash
sbatch run_phase4_single.sbatch
```

Experiments were performed for:

- 10 species
- 20 species
- 30 species
- 50 species
- 55 species

using multiple MPI process counts.

### Multi-node

Create the output directories:

```bash
mkdir -p logs_phase4_nodes results_phase4_nodes
```

Experiments were performed with **32 MPI processes** distributed over different numbers of nodes:

```bash
sbatch --nodes=1 --ntasks=32 --ntasks-per-node=32 run_phase4_nodes.sbatch

sbatch --nodes=2 --ntasks=32 --ntasks-per-node=16 run_phase4_nodes.sbatch

sbatch --nodes=4 --ntasks=32 --ntasks-per-node=8 run_phase4_nodes.sbatch

sbatch --nodes=8 --ntasks=32 --ntasks-per-node=4 run_phase4_nodes.sbatch
```

---

## Output

Each execution generates:

- a distance matrix (`DMat_*`);
- execution logs;
- summary files (`resumen_*.tsv`) containing execution status, Phase 3 time, Phase 4 time, total program time and wall-clock time.
---

## Repository Structure

```text
.
├── include/
├── src/
├── main.cpp
├── Makefile
├── filelist_10
├── filelist_20
├── filelist_30
├── filelist_50
├── filelist_55
├── patterns_clean.txt
├── run_phase4_single.sbatch
├── run_phase4_nodes.sbatch
├── README.md
└── COPYING
```
