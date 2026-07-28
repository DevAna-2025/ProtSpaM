# Prot-SpaM — Paralelización MPI (TFM)

Este repositorio contiene el trabajo desarrollado en el Trabajo Fin de Máster
**"Reconstrucción filogenética de secuencias de proteoma completo en paralelo
sobre sistemas de memoria distribuida"** (Máster Interuniversitario en
Computación de Altas Prestaciones, UDC/USC).

Parte de [Prot-SpaM](https://github.com/jschellh/ProtSpaM), una herramienta de
reconstrucción filogenética libre de alineamiento basada en *spaced-word
matches*, cuya versión original emplea paralelismo de memoria compartida
(OpenMP). Este trabajo desarrolla una versión de **memoria distribuida con MPI**,
capaz de ejecutarse sobre varios nodos de un clúster.

## Estructura de ramas

La rama `master` conserva la **versión original de referencia** (Prot-SpaM con
OpenMP). El trabajo del TFM está organizado en las siguientes ramas:

- **`feat/seq`** — versión secuencial de referencia, empleada como base limpia
  del análisis y como patrón de verificación de corrección.
- **`feat/mpi-phase3-a`** — paralelización MPI de la Fase 3 (lectura
  centralizada, Opción A).
- **`feat/mpi-phase3-b`** — paralelización MPI de la Fase 3 (lectura
  distribuida, Opción B).
- **`feat/mpi-phase4-metacache`** — paralelización MPI de la Fase 4, variante
  con comunicación bloqueante (*metacache*).
- **`feat/mpi-phase4-metacache-isend-calcopt`** — variantes con envíos no
  bloqueantes (*isend*) y con precálculo (*isend_opt*).

Cada rama incluye su propio README con instrucciones de compilación y ejecución.

---

_A continuación se conserva el README original de Prot-SpaM._

---

# Prot-SpaM
### Note
Currently this program only supports Linux.
### Compilation
cd into the root directory (containing the 'Makefile') and type:
`````make ```
### Run
````./protspam [options] -l <filelist> ```
### Filelist
The program takes a plain text file containing the relative paths to each input
dataset. To create your 'filelist' simply type:
``` ls -1 path/to/input/* > filelist ```
This will list each file in specified directory, one file per line.
### Options
```
	-h/-?: print this help and exit
	-w <integer>: pattern weight (default 6)
	-d <integer>: number of don't-care positions (default 40)
	-s <integer>: the minimum score of a spaced-word match to be considered homologous (default: 0)
	-m <integer>: number of patterns used (default 5)
	-t <integer>: number of threads (default: omp_get_max_threads() )
	-o <filename>: filename for distance matrix (default: DMat)
	-l <filename>: specify a list of files to read as input
	-z : if option is set, the pattern set used will be stored in patterns.txt"
	-p <filename>: filename of pattern set to load and reuse"
````
### Sequence format:
Sequence must be in FASTA format. All protein sequences of one proteome must be contained in one FASTA file.
Example:
`````
>Protein1
RAKSDLKEASDKE..
>Protein2
ATSDLAGTASDKE..
>Protein3
ARNCQEFGSDSDW..
..
