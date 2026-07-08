# ProtSpaM-MPI - Fase 4

Extension con MPI de **Prot-SpaM**, desarrollada como parte del TFM en Computacion de Altas Prestaciones.

Esta carpeta contiene la version usada para el avance de la fase 4. 

## Proyecto original

Este trabajo se basa en:

> Leimeister, C. A., Schellhorn, J., Schoebel, M., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).  
> Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.

Repositorio original:

https://github.com/jschellh/ProtSpaM

## Cambios principales de esta fase

### Fase 3

- Las especies se distribuyen entre procesos MPI.
- El reparto se balancea usando el tamano de las secuencias.
- Cada proceso calcula los spaced-words de sus especies locales.

### Fase 4

- Los pares de especies se calculan en paralelo.
- Los spaced-words se procesan patron por patron para reducir memoria.
- En multinodo se evita enviar datos remotos a procesos que no los necesitan.
- Se precalculan posiciones don't-care y bloques de spaced-words con la misma clave para reducir trabajo repetido en `calc_matches`.
- El programa imprime tiempos separados para fase 3, fase 4 y tiempo total.

## Requisitos

En FinisTerrae III se usaron los siguientes modulos:

```bash
module load cesga/2025
module load gcc
module load openmpi/5.0.9
```

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

Los proteomas no se incluyen en el repositorio. Antes de ejecutar, debe existir la carpeta `data/` y los archivos referenciados por cada filelist.

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
```

No se uso `filelist_40` en el informe de avance.

## Ejecucion manual

Ejemplo:

```bash
mpirun -np 4 ./bin/Debug/protspam \
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

Salida principal:

```text
logs_phase4_single_64g/resumen_single_64g_<job_id>.tsv
logs_phase4_single_64g/*.log
results_phase4_single_64g/DMat_*
```

### Multinodo

El benchmark multinodo se ejecuto lanzando el mismo script varias veces, cambiando los parametros de SLURM para reflejar el numero de nodos. El script guarda `SLURM_JOB_NUM_NODES`, `SLURM_NTASKS`, `ntasks_per_node` y `SLURM_NODELIST` en el resumen TSV.

Crear carpetas:

```bash
mkdir -p logs_phase4_nodes_allnp_64g results_phase4_nodes_allnp_64g
```

Comandos usados para reflejar distintos numeros de nodos:

```bash
sbatch --nodes=1 --ntasks=32 --ntasks-per-node=32 run_phase4_nodes.sbatch

sbatch --nodes=2 --ntasks=32 --ntasks-per-node=16 run_phase4_nodes.sbatch

sbatch --nodes=4 --ntasks=32 --ntasks-per-node=8 run_phase4_nodes.sbatch

sbatch --nodes=8 --ntasks=32 --ntasks-per-node=4 run_phase4_nodes.sbatch
```

Configuracion usada:

```text
Memoria: 64 GB por nodo
Procesos MPI evaluados dentro del script: 1, 2, 4, 8, 16, 32
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

Salida principal:

```text
logs_phase4_nodes_allnp_64g/resumen_nodes<N>_allnp_64g_<job_id>.tsv
logs_phase4_nodes_allnp_64g/*.log
results_phase4_nodes_allnp_64g/DMat_*
```

Donde `<N>` corresponde al numero de nodos asignado por SLURM.

## Resultados generados

Cada ejecucion genera:

- una matriz de distancias `DMat_*`;
- un log individual;
- una fila en un TSV resumen con estado, codigo de salida, tiempo de fase 3, tiempo de fase 4, tiempo total y tiempo real de ejecucion.

Los TSV finales usados para el informe de avance corresponden a ejecuciones con 10, 20, 30, 50 y 55 especies, 64 GB de memoria por nodo y 5 repeticiones por configuracion.

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
|-- patterns_clean.txt
|-- run_phase4_single.sbatch
|-- run_phase4_nodes.sbatch
|-- README.md
|-- COPYING
```
