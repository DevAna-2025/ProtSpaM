# ProtSpaM-MPI - Fase 4 (variante: metacache)

Extensión con MPI de **Prot-SpaM**, desarrollada como parte del TFM en Computación de Altas Prestaciones.

Esta rama contiene la variante **metacache** de la Fase 4, base sobre la que se construyen el resto de optimizaciones (`isend`, `metacache_calcopt`, `isend_calcopt`).


## Proyecto original

Este trabajo se basa en:

> Leimeister, C. A., Schellhorn, J., Schoebel, M., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).  
> Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.

Repositorio original:

https://github.com/jschellh/ProtSpaM

## Características:

Esta rama parte de **`feat/mpi-phase3-a`** (lectura centralizada: rank 0 lee
todas las especies y las reparte) como base, y anade la paralelización de la
Fase 4 (cálculo de matches y matriz de distancias):

- Cada especie tiene un proceso propietario; los demas procesos reciben sus
  datos solo si los necesitan para algun par pendiente.
  
- La metadata de cada especie remota (cabecera, secuencia y posiciones de
  inicio) se envia **una sola vez por especie**, antes del bucle de patrones,
  y se guarda en una tabla en memoria (`remote_metadata_cache`) para no
  reenviarla en cada uno de los 5 patrones.
  
- Antes de comunicar nada se calcula que procesos necesitan realmente cada
  especie remota, y solo se envia a esos procesos.
  
- Los spaced-words (que si cambian en cada patron) se recalculan y comunican
  patron a patron, descartando los del patron anterior para no acumular
  memoria.
  
- La comunicacion de metadata y de spaced-words usa `MPI_Send`/`MPI_Recv`
  bloqueante. 


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
reciente (ver `filelists/README.md` para el detalle de su construccion).

## Ejecucion manual

Ejemplo:

```bash
mpirun -np 4 ./bin/Debug/protspam_metacache \
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

### Un nodo

Crear carpetas:

```bash
mkdir -p logs_phase4_single_64g results_phase4_single_64g
```

Ejecutar:

```bash
sbatch run_phase4_single.sbatch
```

Configuracion usada:

```text
Nodos: 1
Memoria: 64 GB
Procesos MPI evaluados: 1, 2, 4, 8, 16, 32
Repeticiones: 5
```

Datasets ejecutados en el script:

```bash
run_species_set 10 filelist_10
run_species_set 20 filelist_20
run_species_set 30 filelist_30
run_species_set 50 filelist_50
run_species_set 55 filelist_55
```

Esta variante se ejecuto ademas, con `filelist_64` (conjunto balanceado, 64
especies), evaluando np = 1, 2, 4, 8, 16, 32 y 64 en un unico nodo (ver
tambien la seccion "Multinodo" para la misma variante ejecutada en varios
nodos).

Salida principal:

```text
logs_phase4_single_64g/resumen_single_64g_<job_id>.tsv
logs_phase4_single_64g/*.log
results_phase4_single_64g/DMat_*
```

### Multinodo

El benchmark multinodo del conjunto balanceado (`filelist_64`) se ejecuto con
una densidad fija de 32 procesos por nodo, variando el numero de nodos: 1, 2,
4 y 8 nodos (equivalentes a 32, 64, 128 y 256 procesos totales). El reparto de
procesos por nodo se controla con `mpirun --map-by ppr:32:node`.

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
Memoria: 64 GB por nodo
Repeticiones: 3
```

Dataset ejecutado en el script: `filelist_64` (conjunto balanceado de 64
especies).

Salida principal:

```text
logs_multinodo/resumen_multinodo_<job_id>.tsv
logs_multinodo/*.log
results_multinodo/DMat_*
```


## Estructura principal

```text
.
|-- include/
|-- src/
|-- data/                          (vacio en el repo; ver Preparacion de datos)
|-- main.cpp
|-- Makefile
|-- filelist_10
|-- filelist_20
|-- filelist_30
|-- filelist_50
|-- filelist_55
|-- filelist_64
|-- patterns_clean.txt
|-- run_phase4_single.sbatch
|-- run_benchmark_multinodo.sbatch
|-- logs_phase4_single_64g/         (logs y TSV resumen; sin DMat_*)
|-- logs_multinodo/                 (logs y TSV resumen; sin DMat_*)
|-- README.md
|-- COPYING
```
