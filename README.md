# ProtSpaM-MPI - Fase  3

Extensión con MPI de **Prot-SpaM**, desarrollada como parte de un Trabajo de Fin de Máster en Computación de Altas Prestaciones.

La implementación actual paraleliza la **etapa de generación de palabras espaciadas (Fase 3)**. Las etapas restantes conservan el comportamiento secuencial original y continúan en desarrollo.

---

## Proyecto Original

Este trabajo está basado en:

> Leimeister, C. A., Schellhorn, J., Schoebel, S., Gerth, M., Bleidorn, C., & Morgenstern, B. (2018).
> *Prot-SpaM: Fast alignment-free phylogeny reconstruction based on whole-proteome sequences.*
> bioRxiv, 306142.

Repositorio original:

https://github.com/jschellh/ProtSpaM

---

## Características

* Paralelización de la Fase 3 (generación de palabras espaciadas) mediante MPI.
* Computación distribuida entre múltiples procesos MPI.
* Experimentos reproducibles utilizando archivos de patrones fijos.
* Compatible con el formato de entrada y los conjuntos de datos originales de Prot-SpaM.

---

## Compilación

```bash
make
```

Ejecutable:

```bash
./bin/Debug/protspam
```

---

## Preparación de los datos

Antes de ejecutar el programa, crea un directorio `data/` y coloca los archivos FASTA según las rutas especificadas en el `filelist` correspondiente.

Los conjuntos de datos FASTA no están incluidos en este repositorio. Los archivos listados en `filelist_10`, `filelist_20` y `filelist_30` deben existir antes de ejecutar los experimentos.

Ejemplo:

```bash
mkdir data
cp /ruta/a/proteomas/*.faa data/
```

Si `filelist_10` contiene:

```text
data/species1.faa
data/species2.faa
data/species3.faa
```

entonces los archivos deben existir en esas ubicaciones.

---

## Ejecución

```bash
mpirun -np <procesos> ./bin/Debug/protspam [opciones] -l <filelist> -p <patrones>
```

Ejemplo:

```bash
mpirun -np 32 ./bin/Debug/protspam \
    -l filelist_30 \
    -p patterns_clean.txt \
    -o DMat_30sp
```

---

## Script de Benchmark

Los experimentos de la **Fase 3 - Opción A** pueden ejecutarse mediante el script de SLURM:

```bash
mkdir -p logs results
sbatch run_opcionA.sbatch
```

El script ejecuta 5 repeticiones para 10, 20 y 30 especies utilizando 1, 2, 4, 8, 16 y 32 procesos MPI. Los registros individuales de ejecución se almacenan en `logs/`, las matrices de distancia en `results/` y se genera un archivo resumen en formato TSV con el estado, el código de salida y el tiempo de ejecución de cada experimento.

Los experimentos están diseñados para ejecutarse dentro de un único nodo. El script solicita 32 tareas MPI y utiliza hasta 32 procesos por ejecución.

---

## Filelists

Los experimentos de este repositorio utilizan:

* `filelist_10`
* `filelist_20`
* `filelist_30`

---

## Patrones

Los experimentos utilizan el archivo de patrones fijo:

```text
patterns_clean.txt
```

con los parámetros por defecto de Prot-SpaM:

* Peso: 6
* Posiciones *don't-care*: 40
* Umbral: 0
* Número de patrones: 5

---

## Estructura del repositorio

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
