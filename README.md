# ProtSpaM-MPI

Extensión con MPI de **Prot-SpaM**, desarrollada como parte de un Trabajo de Fin de Máster en Computación de Altas Prestaciones.

Esta implementación paraleliza las etapas más costosas de Prot-SpaM, al tiempo que reduce el consumo de memoria durante el cálculo de coincidencias mediante una tubería (*pipeline*) de procesamiento a nivel de patrón.

## Proyecto Original

Basado en:

> Leimeister, C. A., Schellhorn, J., Schoebel, S., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).  
> *Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.*

Repositorio original:

https://github.com/jschellh/ProtSpaM

---

## Implementación MPI

### Fase 3 – Generación de palabras espaciadas

- Las especies se distribuyen entre los distintos *ranks* de MPI.
- El balanceo de carga se realiza en función de la longitud de las secuencias.
- Cada *rank* calcula las palabras espaciadas correspondientes a las especies que tiene asignadas.

### Fase 4 – Cálculo de coincidencias

- Los pares de especies se procesan en paralelo.
- Las palabras espaciadas se transmiten (*streaming*) un patrón cada vez para reducir el uso de memoria.
- Las palabras espaciadas remotas se intercambian bajo demanda, en lugar de replicar el conjunto de datos completo.
- Las longitudes de secuencias consecutivas con la misma clave (*equal-key run lengths*) y las posiciones *don't-care* se precalculan para acelerar el proceso de búsqueda de coincidencias.
---


## Compilación

```bash
make clean
make
```

Ejecutable:

```text
./bin/Debug/protspam
```

---

## Preparación de los datos

Crea el directorio de datos y copia los archivos de proteomas referenciados por los *filelists*.

```bash
mkdir -p data
cp /ruta/a/proteomas/*.faa data/
```

Los conjuntos de datos **no están incluidos** en este repositorio.

Los experimentos utilizaron los siguientes *filelists*:

- `filelist_10`
- `filelist_20`
- `filelist_30`
- `filelist_50`
- `filelist_55`

---

## Ejecución

Ejemplo:

```bash
mpirun -np 4 ./bin/Debug/protspam \
    -l filelist_20 \
    -p patterns_clean.txt \
    -o DMat_20sp_np4
```

El programa muestra:

```text
Tiempo de generación de palabras espaciadas
Tiempo de cálculo de coincidencias
Tiempo total
```

---

## Benchmarks

### Un solo nodo

Crea los directorios de salida:

```bash
mkdir -p logs_phase4_single_64g results_phase4_single_64g
```

Ejecuta:

```bash
sbatch run_phase4_single.sbatch
```

Los experimentos se realizaron para:

- 10 especies
- 20 especies
- 30 especies
- 50 especies
- 55 especies

utilizando diferentes cantidades de procesos MPI.

### Múltiples nodos

Crea los directorios de salida:

```bash
mkdir -p logs_phase4_nodes results_phase4_nodes
```

Los experimentos se realizaron con **32 procesos MPI** distribuidos en distintos números de nodos:

```bash
sbatch --nodes=1 --ntasks=32 --ntasks-per-node=32 run_phase4_nodes.sbatch

sbatch --nodes=2 --ntasks=32 --ntasks-per-node=16 run_phase4_nodes.sbatch

sbatch --nodes=4 --ntasks=32 --ntasks-per-node=8 run_phase4_nodes.sbatch

sbatch --nodes=8 --ntasks=32 --ntasks-per-node=4 run_phase4_nodes.sbatch
```

---

## Salida

Cada ejecución genera:

- una matriz de distancias (`DMat_*`);
- registros de ejecución (*logs*);
- archivos resumen (`resumen_*.tsv`) que contienen el estado de la ejecución, el tiempo de la Fase 3, el tiempo de la Fase 4, el tiempo total del programa y el tiempo de ejecución real (*wall-clock time*).

---

## Estructura del repositorio

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
