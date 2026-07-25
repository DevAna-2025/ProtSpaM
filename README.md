# ProtSpaM-MPI — Fase 4 (variante isend: comunicación no bloqueante)

Extensión con MPI de **Prot-SpaM**, desarrollada como parte de un Trabajo de
Fin de Máster en Computación de Altas Prestaciones.

Esta rama paraleliza mediante MPI la **Fase 4** (cálculo de las coincidencias y
de la matriz de distancias). Mantiene exactamente el mismo diseño que la
variante *metacache* —reparto, matriz de necesidad remota, caché de metadatos y
streaming por patrón— y sustituye únicamente el mecanismo de envío: en lugar de
un `MPI_Send` bloqueante por destino, lanza todos los `MPI_Isend` de una tanda
sin esperar entre ellos y espera una sola vez con `MPI_Waitall`. Así, el proceso
propietario no queda serializado esperando turno destino a destino.

## Cambios respecto a Prot-SpaM original

- Reparto estático de las especies por bloques entre los procesos MPI.
- Cálculo de la matriz de distancias distribuido: cada proceso calcula los
  pares en los que interviene alguna de sus especies locales, comunicando las
  palabras espaciadas de las especies remotas que necesita.
- **Matriz de necesidad remota** (`MPI_Allgather`): el propietario de una
  especie solo la envía a los procesos que realmente la van a usar.
- **Caché de metadatos remotos**: el encabezado, la secuencia y los metadatos
  de cada especie remota se comunican una sola vez, antes del bucle de patrones.
- **Streaming patrón a patrón**: las palabras espaciadas se calculan y comunican
  patrón a patrón, reduciendo el pico de memoria en un factor igual al número de
  patrones.
- **Comunicación no bloqueante** (`MPI_Isend` + `MPI_Waitall`): los envíos de una
  tanda se lanzan sin esperar entre sí y se espera una sola vez al final. Los
  buffers de cada tanda (`PendingSend`) se mantienen vivos hasta completar el
  `MPI_Waitall` correspondiente.
- Carga de patrones desde fichero fijo para garantizar la reproducibilidad.

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
