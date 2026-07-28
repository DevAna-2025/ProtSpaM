# ProtSpaM (Versión Secuencial en C++)

Esta es una adaptación **secuencial** de
[ProtSpaM](https://github.com/jschellh/ProtSpaM), una herramienta para estimar
distancias filogenéticas entre proteínas basada en **spaced-word matches**.

La versión original incluye paralelismo con **OpenMP**. Esta versión elimina
toda la paralelización para ejecutarse de manera **estrictamente secuencial**,
manteniendo la misma funcionalidad principal. Sirve como **implementación de
referencia**: una base limpia para el análisis del algoritmo y el patrón contra
el que se verifica que las versiones paralelas producen el mismo resultado.

---

## Cambios respecto a la versión original

- Se creó un archivo `main.cpp`  que es la versión secuencial del original
- Se eliminaron las dependencias de **OpenMP** (`#include <omp.h>`,
  `omp_get_wtime`, `omp_set_num_threads`, `#pragma omp parallel for`).
- El sistema de medición de tiempo ahora utiliza `std::chrono` en lugar de
  funciones de OpenMP.
- Se eliminó el parámetro `-t` (número de hilos), innecesario en código
  secuencial.
- Se añadió la rutina `write_words_checksum`, que vuelca a fichero las palabras
  espaciadas generadas. Sirve como referencia para verificar que las versiones
  paralelas producen exactamente el mismo resultado.

---

## Compilación

1. Necesitas una carpeta `data/` con todos los archivos FASTA.
2. Revisa que las rutas dentro del archivo `filelist` existan en `data/`.
3. Compila con el `Makefile`:

```bash
make
```

Genera el ejecutable en `./bin/Debug/protspam`.

---

## Uso

El programa acepta los mismos parámetros que la versión original. Ejemplo con
múltiples archivos de entrada:

```bash
./bin/Debug/protspam -w 6 -d 40 -m 5 -l filelist -p patterns.txt
```

### Opciones principales

- `-w <int>` : Peso del patrón (**default:** 6)
- `-d <int>` : Número de posiciones "don't-care" (**default:** 40)
- `-s <int>` : Umbral mínimo para considerar un spaced-word match como homólogo (**default:** 0)
- `-m <int>` : Número de patrones (**default:** 5)
- `-o <file>` : Nombre del archivo de salida con la matriz de distancias (**default:** DMat)
- `-l <file>` : Lista de archivos de entrada en formato multifasta
- `-p <file>` : Cargar conjunto de patrones predefinidos
- `-z` : Guardar patrones generados en `patterns.txt`
- `-r` : Generar puntajes individuales para cada par de secuencias (spamogramas)

---

## Salida

El programa genera:

- **Matriz de distancias** → archivo de salida (`DMat` por defecto).
- **Spamogramas (scores)** → en el directorio `scores/` si se usa la opción `-r`.

---

## Contacto

Herramienta original: jendrik.schellhorn@stud.uni-goettingen.de
