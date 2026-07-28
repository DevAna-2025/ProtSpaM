# ProtSpaM-MPI — Fase 4 (variante isend_opt: precálculo por patrón)

Extensión con MPI de **Prot-SpaM**, desarrollada como parte de un Trabajo de
Fin de Máster en Computación de Altas Prestaciones.

Esta rama parte de la variante *isend* y añade una optimización sobre la
**ruta de referencia secuencial** (`np = 1`): traslada a la función
`calc_matches` el precálculo por patrón de las posiciones *don't-care* y de los
bloques de palabras espaciadas con la misma clave. La ruta paralela ya
incorporaba ese precálculo desde la variante *metacache*, por lo que el esquema
de comunicación no cambia respecto a *isend*.

El efecto del cambio solo puede observarse con `np = 1`, ya que con más
procesos el programa ejecuta la ruta paralela. La optimización se evaluó con la
hipótesis de reducir trabajo repetido en la referencia secuencial. Sin embargo,
en los ensayos recogidos en el TFM no produjo una mejora: para `np = 1`,
`isend_opt` registró 3755,3 s frente a 3141,4 s de `isend`. Se conserva como
resultado experimental y como referencia para futuros ajustes.

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
  patrón a patrón, reduciendo el pico de memoria.
- **Comunicación no bloqueante** (`MPI_Isend` + `MPI_Waitall`): los envíos de
  una tanda se publican antes de realizar una espera conjunta. Los búferes se
  mantienen vivos hasta completar `MPI_Waitall`. Este mecanismo no garantiza,
  por sí solo, solapamiento entre comunicación y cómputo.
- **Precálculo por patrón en `calc_matches`**: las posiciones *don't-care* y
  los bloques de palabras con la misma clave se calculan una sola vez por
  patrón, en lugar de recalcularse en cada comparación. Solo afecta a la ruta
  secuencial (`np = 1`).
- Carga de patrones desde fichero fijo para garantizar la reproducibilidad.

## Proyecto original

> Leimeister, C. A., Schellhorn, J., Dörrer, S., Gerth, M., Bleidorn, C.,
> & Morgenstern, B. (2019). *Prot-SpaM: fast alignment-free phylogeny
> reconstruction based on whole-proteome sequences.* GigaScience, 8(3), giy148.

Repositorio original: https://github.com/jschellh/ProtSpaM

## Compilación

```bash
make
