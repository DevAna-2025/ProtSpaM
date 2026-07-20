# ProtSpaM-MPI - Fase 4 (variante: isend_calcopt)

Extensión con MPI de **Prot-SpaM**, desarrollada como parte del TFM en Computación de Altas Prestaciones.

Esta rama contiene la variante **isend_calcopt** de la Fase 4: parte de `isend` y optimiza ademas el calculo de matches (`calc_matches`) usado en la ejecucion con un unico proceso.


## Proyecto original

Este trabajo se basa en:

> Leimeister, C. A., Schellhorn, J., Schoebel, M., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).  
> Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.

Repositorio original:

https://github.com/jschellh/ProtSpaM

## Características:

Esta rama parte de **`isend`** (rama `feat/mpi-phase4-metacache-isend`) como base, que
a su vez parte de `metacache` (`feat/mpi-phase4-metacache`), que a su vez
parte de `feat/mpi-phase3-a` (lectura centralizada). De `isend` se conserva
todo: la metadata en cache por especie, el envio selectivo, el streaming de
spaced-words por patron y el envio no bloqueante (`MPI_Isend` + `MPI_Waitall`)
al comunicar datos a varios destinos (ver los README de `metacache` e `isend`
para el detalle).

El cambio en esta rama esta en `calc_matches` (`calc_matches.cpp`), la funcion
que compara un par de especies completo. Esta funcion solo se usa en la
ejecucion con un unico proceso (`size == 1`); con mas de un proceso, el
pipeline MPI ya usaba una logica equivalente de forma interna. Las
optimizaciones que incorpora `calc_matches` en esta rama son:

- **Posiciones don't-care precalculadas** (`dontcare_positions`): en el
  patron original, todas las posiciones se recorren para calcular el score;
  aqui se calculan una vez, al principio de cada patron, las posiciones
  marcadas como don't-care, que son las unicas relevantes para el calculo de
  matches.

- **Comparacion restringida a esas posiciones** (`score_dontcare_positions`):
  el score y los mismatches entre dos spaced-words se calculan recorriendo
  solo esas posiciones precalculadas, en lugar de todo el patron.

- **Longitud de bloques con la misma key precalculada** (`key_run_lengths`):
  los spaced-words estan ordenados por `key`; esta tabla indica, para cada
  posicion, cuantos elementos consecutivos comparten la misma key. Mantiene
  el mismo comportamiento que la funcion `multiMatch` original, pero evita
  recalcular esa busqueda repetidamente dentro de los bucles principales.

- El resultado del calculo (tasa de mismatches y matriz de distancias) es el
  mismo que en `isend`; el cambio reduce trabajo repetido dentro de
  `calc_matches`, no altera la logica de comparacion.


## Compilacion

```bash
make clean
make
```

El ejecutable queda en:

```text
./bin/Debug/protspam
```

## Preparacion de datos

Los conjuntos de datos FASTA no están incluidos en este repositorio, fueron utilizados los archivos referenciados en el repositorio original, veáse el enlace del conjunto de datos [aquí](http://projects.gobics.de/data/protspam/paperData.tgz) . Antes de ejecutar, debe existir la carpeta `data/` y los archivos referenciados por cada filelist.

Ejemplo:

```bash
mkdir -p data
cp /ruta/a/proteomas/*.faa data/
```

Los experimentos de fase 4 se prepararon con estos filelists:

```text
filelist_10
filelist_20
filelist_30
filelist_50
filelist_55
filelist_64
```

No se uso `filelist_40` en el informe de avance. `filelist_64` es un conjunto
balanceado de 64 especies (sin la especie de mayor tamaño ni las de menor
tamano del conjunto original), usado en el benchmark de escalabilidad mas
reciente ( [ver](https://github.com/DevAna-2025/ProtSpaM/blob/feat/mpi-phase4-metacache/README_FILELIST.md))  para el detalle de su construccion).

## Ejecucion manual

Ejemplo:

```bash
mpirun -np 4 ./bin/Debug/protspam_block_pipeline_metacache_isend_calcopt \
    -l filelist_20 \
    -p patterns_clean.txt \
    -o DMat_20sp_np4
```

La salida del programa incluye:

```text
Tiempo spaced-words
Tiempo matches
Tiempo total
```

## Benchmarks

Antes de ejecutar los batch scripts es necesario crear las carpetas donde se escribiran logs y matrices de salida. Si las carpetas no existen, SLURM puede fallar al abrir los archivos indicados en `#SBATCH --output` y `#SBATCH --error`.

### Un nodo (filelist_55)

Al igual que `isend`, esta variante no se ejecuto con el script generico
`run_variant_55.sbatch`, sino con un script propio
(`run_isend_calcopt_55.sbatch`). Este script compara especificamente esta
variante frente a `isend` (version base) para evaluar si la optimizacion de
`calc_matches` aporta mejora adicional sobre el fix de comunicacion. Solo
evalua np = 16 y 32 (no 1, 2, 4, 8, 64), con 3 repeticiones.

Crear carpetas (las crea tambien el propio script si no existen):

```bash
mkdir -p logs_isend_55 results_isend_55
```

Ejecutar:

```bash
sbatch run_isend_calcopt_55.sbatch
```

Configuracion usada:

```text
Nodos: 1
Memoria: 64 GB
Filelist: filelist_55
Procesos MPI evaluados: 16, 32
Repeticiones: 3
```

Salida principal:

```text
logs_isend_55/resumen_isend_calcopt_<job_id>.tsv
logs_isend_55/isend_calcopt_55sp_np*_rep*_<job_id>.log
results_isend_55/DMat_isend_calcopt_55sp_np*_rep*_<job_id>
```

Nota: esta rama comparte las carpetas `logs_isend_55/` y `results_isend_55/`
con la rama `isend` (ambos scripts escriben en las mismas carpetas); los
ficheros se distinguen por el prefijo `isend_calcopt_` frente a
`isend_metacache_`.

### Un nodo (filelist_64)

Esta variante se ejecuto ademas, con `filelist_64` (conjunto balanceado, 64
especies), evaluando np = 1, 2, 4, 8, 16, 32 y 64 en un unico nodo (ver
tambien la seccion "Multinodo" para la misma variante ejecutada en varios
nodos).

El script `run_benchmark_64.sbatch` evalua en una misma ejecucion las
variantes `metacache`, `isend` e `isend_opt` (esta rama, misma optimizacion
de `calc_matches` descrita arriba, identificada como `isend_opt` en dicho
script), para optimizar el uso del cluster.El TSV resumen incluidos en esta rama corresponden unicamente a la variante
`isend_opt` (ficheros con prefijo `isend_opt_`, columna `variante = isend_opt`
en el TSV); las otras variantes se documentan en sus ramas correspondientes.

Crear carpetas:

```bash
mkdir -p logs_bench64 results_bench64
```

Ejecutar:

```bash
sbatch run_benchmark_64.sbatch
```

Configuracion usada:

```text
Nodos: 1
Memoria: 64 GB
Filelist: filelist_64
Procesos MPI evaluados: 1, 2, 4, 8, 16, 32, 64
Repeticiones: 3
```

Salida principal (filtrada a esta variante, prefijo `isend_opt_`):

```text
logs_bench64/resumen_bench64_<job_id>.tsv
logs_bench64/isend_opt_64sp_np*_rep*_<job_id>.log
results_bench64/DMat_isend_opt_64sp_np*_rep*_<job_id>
```

### Multinodo

El benchmark multinodo del conjunto balanceado (`filelist_64`) se ejecuto con
una densidad fija de 32 procesos por nodo, variando el numero de nodos: 1, 2,
4 y 8 nodos (equivalentes a 32, 64, 128 y 256 procesos totales). El reparto de
procesos por nodo se controla con `mpirun --map-by ppr:32:node`.

El script `run_benchmark_multinodo.sbatch` evalua en una misma ejecucion las
variantes `metacache`, `isend` e `isend_opt` (esta rama, identificada como
`isend_opt` en dicho script), para optimizar el uso del cluster.El TSV resumen incluidos en esta rama corresponden
unicamente a la variante `isend_opt` (ficheros con prefijo `isend_opt_`,
columna `variante = isend_opt` en el TSV); las otras variantes se documentan
en sus ramas correspondientes.

Crear carpetas:

```bash
mkdir -p logs_multinodo results_multinodo
```

Ejecutar:

```bash
sbatch run_benchmark_multinodo.sbatch
```

Configuracion usada:

```text
Nodos evaluados: 1, 2, 4, 8
Procesos por nodo (PPN): 32
Procesos totales evaluados: 32, 64, 128, 256
Memoria: 64 GB por nodo
Filelist: filelist_64
Repeticiones: 3
```

Salida principal (filtrada a esta variante, prefijo `isend_opt_`):

```text
logs_multinodo/resumen_multinodo_<job_id>.tsv
logs_multinodo/isend_opt_64sp_*nodos_np*_rep*_<job_id>.log
results_multinodo/DMat_isend_opt_64sp_*nodos_np*_rep*_<job_id>
```


## Estructura principal

```text
.
|-- include/
|-- src/
|-- data/                        
|-- main.cpp
|-- Makefile
|-- filelist_10
|-- filelist_20
|-- filelist_30
|-- filelist_50
|-- filelist_55
|-- filelist_64
|-- patterns_clean.txt
|-- run_benchmark_64.sbatch
|-- run_benchmark_multinodo.sbatch
|-- README.md
|-- README_FILELIST.md
|-- COPYING
```
