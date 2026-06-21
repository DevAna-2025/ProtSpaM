# ProtSpaM-MPI

This repository contains an MPI extension of **Prot-SpaM** developed as part of a Master's Thesis in High Performance Computing.

The current implementation focuses on the parallelization of the **spaced-word generation stage (Phase 3)** of the original workflow. The remaining stages preserve the original sequential behaviour and are under development.

---

## Original Project

This work is based on:

> Leimeister, C. A., Schellhorn, J., Schoebel, S., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).
> *Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.*
> bioRxiv, 306142.

Original repository:

https://github.com/jschellh/ProtSpaM

---

## Current Features

* MPI-based parallelization of the spaced-word generation stage (Phase 3).
* Distributed computation of local spaced-words across multiple MPI processes.
* Reproducible experiments using fixed pattern files.
* Compatibility with the original Prot-SpaM input format and datasets.

---

## Work in Progress

* Parallelization of the matches computation stage (Phase 4).
* Reduction of unnecessary data movement between phases.

---

## Compilation

From the root directory:

```bash
make
```

The executable will be generated as:

```bash
./bin/Debug/protspam
```

---

## Running the Program

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

## Input Datasets

The program expects protein datasets in multi-FASTA format, where each FASTA file contains all proteins belonging to a single species/proteome.

Datasets used in the original Prot-SpaM publication can be downloaded from:

http://projects.gobics.de/data/protspam/paperData.tgz

---

## Filelist

Input files are specified through a plain text file containing one FASTA file per line.

To automatically create a filelist:

```bash
ls -1 path/to/input/* > filelist
```

For the experiments presented in this repository:

* `filelist_10`
* `filelist_20`
* `filelist_30`

---

## Pattern Files

Prot-SpaM can either generate patterns randomly or load them from a file.

For reproducible experiments, this project uses:

```text
patterns_clean.txt
```

with the default Prot-SpaM parameters:

* Weight: 6
* Number of don't-care positions: 40
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



Partial MPI parallelization of Prot-SpaM developed as part of a Master's Thesis in High Performance Computing.
