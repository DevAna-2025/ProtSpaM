# ProtSpaM-MPI

MPI extension of **Prot-SpaM** developed as part of a Master's Thesis in High Performance Computing.

The current implementation parallelizes the **spaced-word generation stage (Phase 3)**. The remaining stages preserve the original sequential behaviour and are under development.

---

## Original Project

This work is based on:

> Leimeister, C. A., Schellhorn, J., Schoebel, S., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).
> *Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.*
> bioRxiv, 306142.

Original repository:

https://github.com/jschellh/ProtSpaM

---

## Features

* MPI-based parallelization of Phase 3 (spaced-words generation).
* Distributed computation across multiple MPI processes.
* Reproducible experiments using fixed pattern files.
* Compatible with the original Prot-SpaM input format and datasets.

---

## Compilation

```bash
make
```

Executable:

```bash
./bin/Debug/protspam
```

---

## Preparing the Data

Before running the program, create a `data/` directory and place the FASTA files according to the paths specified in the corresponding `filelist`.

The FASTA datasets are not included in this repository. The files listed in `filelist_10`, `filelist_20` and `filelist_30` must exist before launching the experiments.

Example:

```bash
mkdir data
cp /path/to/proteomes/*.faa data/
```

If `filelist_10` contains:

```text
data/species1.faa
data/species2.faa
data/species3.faa
```

then the files must exist at those locations.

---

## Running

```bash
mpirun -np <processes> ./bin/Debug/protspam [options] -l <filelist> -p <patterns>
```

Example:

```bash
mpirun -np 32 ./bin/Debug/protspam \
    -l filelist_30 \
    -p patterns_clean.txt \
    -o DMat_30sp
```

---

## Benchmark Script

The Phase 3 Option A experiments can be launched with the SLURM script:

```bash
mkdir -p logs results
sbatch run_opcionA.sbatch
```

The script runs 5 repetitions for 10, 20 and 30 species using 1, 2, 4, 8, 16 and 32 MPI processes. It writes individual execution logs to `logs/`, distance matrices to `results/`, and a summary TSV file with the status, exit code and elapsed time of each run.

The experiments are intended to be executed within a single node. The script requests 32 MPI tasks and uses up to 32 processes per run.

---

## Filelists

The experiments in this repository use:

* `filelist_10`
* `filelist_20`
* `filelist_30`

---

## Patterns

Experiments use the fixed pattern file:

```text
patterns_clean.txt
```

with the default Prot-SpaM parameters:

* Weight: 6
* Don't-care positions: 40
* Threshold: 0
* Number of patterns: 5

---

## Repository Structure

```text
.
|-- data/
|-- filelist_10
|-- filelist_20
|-- filelist_30
|-- patterns_clean.txt
|-- main.cpp
|-- run_opcionA.sbatch
|-- logs/
|-- results/
|-- src/
|-- include/
|-- Makefile
|-- README.md
```
