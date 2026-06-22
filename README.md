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
mpirun -np <processes> ./bin/Debug/protspam  [options] -l <filelist> -p <patterns>
```

Example:

```bash
mpirun -np 32 ./bin/Debug/protspam  \
    -l filelist_30 \
    -p patterns_clean.txt \
    -o DMat_30sp
```

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
├── data/
├── filelist_10
├── filelist_20
├── filelist_30
├── patterns_clean.txt
├── main.cpp
├── logs/
├── results/
├── src/
├── Makefile
└── README.md
```
