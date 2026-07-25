# ProtSpaM-MPI — Fase 3 (Opción A: lectura centralizada)

Extensión con MPI de **Prot-SpaM**, desarrollada como parte de un Trabajo de
Fin de Máster en Computación de Altas Prestaciones.

Esta rama paraleliza mediante MPI la **Fase 3** (generación de las palabras
espaciadas). En esta variante, el proceso 0 lee todas las secuencias del disco
y reparte a cada proceso su bloque de especies mediante comunicación punto a
punto. Las fases restantes conservan el comportamiento secuencial original.

## Cambios respecto a Prot-SpaM original

- Reparto estático de las especies por bloques entre los procesos MPI.
- Cálculo de las palabras espaciadas distribuido: cada proceso opera sobre sus
  especies locales, sin comunicación durante el cómputo.
- **Lectura centralizada**: solo el proceso 0 accede al disco y distribuye las
  secuencias al resto (`send_species` / `recv_species`).
- Carga de patrones desde fichero fijo para garantizar la reproducibilidad
  (se prescinde de la generación probabilística por defecto).

## Proyecto original

> Leimeister, C. A., Schellhorn, J., Dörrer, S., Gerth, M., Bleidorn, C.,
> & Morgenstern, B. (2019). *Prot-SpaM: fast alignment-free phylogeny
> reconstruction based on whole-proteome sequences.* GigaScience, 8(3), giy148.

Repositorio original: https://github.com/jschellh/ProtSpaM

## Compilación

```bash
make
```

Genera el ejecutable en `./bin/Debug/protspam`.

## Ejecución

```bash
mpirun -np <procesos> ./bin/Debug/protspam -l <filelist> -p <patrones> -o <salida>
```

Ejemplo:

```bash
mpirun -np 32 ./bin/Debug/protspam -l filelist -p patterns_clean.txt -o DMat
```

Donde:

- `<filelist>`: fichero de texto con la ruta a cada FASTA, una por línea.
- `<patrones>`: fichero de patrones fijos.
- `<salida>`: fichero de la matriz de distancias resultante (formato PHYLIP).

## Parámetros

Los patrones se cargan desde `patterns_clean.txt`, con la configuración por
defecto de Prot-SpaM:

- Peso del patrón: 6
- Posiciones *don't-care*: 40
- Umbral: 0
- Número de patrones: 5

## Datos de entrada

Los ficheros FASTA no se incluyen en el repositorio. Cada especie es un fichero
FASTA con su proteoma completo, y el `filelist` contiene la ruta a cada uno,
una por línea:

```text
data/species1.faa
data/species2.faa
data/species3.faa
```
Se emplearon los ficheros referenciados en el repositorio original de
Prot-SpaM, disponibles en
http://projects.gobics.de/data/protspam/paperData.tgz

Todas las rutas listadas deben existir antes de ejecutar el programa.
