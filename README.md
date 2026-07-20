# ProtSpaM-MPI - Fase 4 (variante: isend)

Extensión con MPI de **Prot-SpaM**, desarrollada como parte del TFM en Computación de Altas Prestaciones.

Esta rama contiene la variante **isend** de la Fase 4: parte de `metacache` y sustituye su comunicación bloqueante por un envío no bloqueante.


## Proyecto original

Este trabajo se basa en:

> Leimeister, C. A., Schellhorn, J., Schoebel, M., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).  
> Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.

Repositorio original:

https://github.com/jschellh/ProtSpaM

## Características:

Esta rama parte de **`metacache`** (rama `feat/mpi-phase4-metacache`) como
base, que a su vez parte de `feat/mpi-phase3-a` (lectura centralizada). De
metacache se conserva todo: la metadata en cache por especie, el envio
selectivo y el streaming de spaced-words por patron (ver el README de
`metacache` para el detalle). El unico cambio es como se envian esos datos:

- **Problema en metacache:** cuando un proceso propietario de una especie
  tenia que enviarla a varios destinos, usaba `MPI_Send` bloqueante uno detras
  de otro. Con mensajes grandes, MPI espera a que el receptor este listo antes
  de devolver el control (rendezvous), asi que el propietario quedaba
  esperando a cada destino en fila antes de pasar al siguiente.

- **Cambio en isend:** los envios a todos los destinos de una misma tanda se
  lanzan con `MPI_Isend` (no bloqueante), sin esperar entre ellos, y se hace
  un unico `MPI_Waitall` al final para confirmar que todos terminaron. Asi el
  progreso de los envios ocurre en paralelo en lugar de en fila.

- Este cambio se aplica en los dos puntos donde metacache enviaba datos a
  varios destinos: el envio de metadata de una especie remota y el envio de
  sus spaced-words por patron.

- El resto del programa (reparto de especies, calculo de matches, fases y
  tiempos reportados) es identico a `metacache`.


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
mpirun -np 4 ./bin/Debug/protspam_block_pipeline_metacache_isend_calcbase \
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

A diferencia de las demas variantes, `isend` no se ejecuto con el script
generico `run_variant_55.sbatch`, sino con un script propio
(`run_isend_55_64g.sbatch`) para ahorrar tiempo de cluster. Este script solo
evalua np = 16, 32 y 64 (no 1, 2, 4, 8), ya que el objetivo de esta corrida
era comparar directamente contra `metacache` en los puntos donde el
desbalance de `filelist_55` es mas relevante.

Crear carpetas (las crea tambien el propio script si no existen):

```bash
mkdir -p logs_isend_55 results_isend_55
```

Ejecutar:

```bash
sbatch run_isend_55_64g.sbatch
```

Configuracion usada:

```text
Nodos: 1
Memoria: 64 GB
Filelist: filelist_55
Procesos MPI evaluados: 16, 32, 64
Repeticiones: 3
```

Salida principal:

```text
logs_isend_55/resumen_isend_<job_id>.tsv
logs_isend_55/isend_metacache_55sp_np*_rep*_<job_id>.log
results_isend_55/DMat_isend_metacache_55sp_np*_rep*_<job_id>
```

### Un nodo (filelist_64)

Esta variante se ejecuto ademas, con `filelist_64` (conjunto balanceado, 64
especies), evaluando np = 1, 2, 4, 8, 16, 32 y 64 en un unico nodo (ver
tambien la seccion "Multinodo" para la misma variante ejecutada en varios
nodos).

El script `run_benchmark_64.sbatch` evalua en una misma ejecucion las
variantes `metacache`, `isend` e `isend_opt`, para optimizar el uso del
cluster. 

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

Salida principal (filtrada a esta variante, prefijo `isend_`):

```text
logs_bench64/resumen_bench64_<job_id>.tsv
logs_bench64/isend_64sp_np*_rep*_<job_id>.log
results_bench64/DMat_isend_64sp_np*_rep*_<job_id>
```

### Multinodo

El benchmark multinodo del conjunto balanceado (`filelist_64`) se ejecuto con
una densidad fija de 32 procesos por nodo, variando el numero de nodos: 1, 2,
4 y 8 nodos (equivalentes a 32, 64, 128 y 256 procesos totales). El reparto de
procesos por nodo se controla con `mpirun --map-by ppr:32:node`.

El script `run_benchmark_multinodo.sbatch` evalua en una misma ejecucion las
variantes `metacache`, `isend` e `isend_opt`, para optimizar el uso del
cluster. 

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

Salida principal (filtrada a esta variante, prefijo `isend_`):

```text
logs_multinodo/resumen_multinodo_<job_id>.tsv
logs_multinodo/isend_64sp_*nodos_np*_rep*_<job_id>.log
results_multinodo/DMat_isend_64sp_*nodos_np*_rep*_<job_id>
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
|-- run_isend_55_64g.sbatch
|-- run_benchmark_64.sbatch
|-- run_benchmark_multinodo.sbatch
|-- README.md
|-- README_FILELIST.md
|-- COPYING
```
