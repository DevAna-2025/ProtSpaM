# Filelists

Los filelists son listas de texto plano con las rutas a los proteomas usados
en cada experimento. Definen que especies se incluyen y en que orden, pero
no contienen las secuencias en si.

## Origen de los proteomas

Los proteomas usados en los filelists `filelist_10`, `filelist_20`,
`filelist_30`, `filelist_50` y `filelist_55` provienen del conjunto de datos
del articulo original de ProtSpaM:

http://projects.gobics.de/data/protspam/paperData.tgz

Los proteomas no se incluyen en este repositorio (ver seccion "Preparacion
de datos" en el README de cada rama). Para reproducir los experimentos,
descargar el `.tgz` anterior y extraerlo en la carpeta `data/`.

## filelist_10 / filelist_20 / filelist_30

Subconjuntos anidados y crecientes de especies (10, 20 y 30), usados en la
evaluacion de la Fase 3 (ramas `feat/mpi-phase3-a` y `feat/mpi-phase3-b`)
para caracterizar la escalabilidad al aumentar progresivamente el tamano del
problema. Todos incluyen `Homo.faa`.

## filelist_50 / filelist_55

Conjuntos usados en la evaluacion de la Fase 4 (variantes `metacache`,
`isend`, `metacache_calcopt`, `isend_calcopt`). `filelist_55` incluye
`Homo.faa`, la especie de mayor tamano del conjunto (75,8 millones de
aminoacidos), junto con las especies de menor tamano (proteomas de
*Wolbachia*, ~0,3 millones de aminoacidos). Esta diferencia de tamano entre
especies (ratio maximo/minimo ≈ 220x) provoca un desbalance de carga entre
procesos, analizado en la memoria del TFM.

## filelist_64

Conjunto balanceado de 64 especies, construido a partir de `filelist_55`
para evaluar la escalabilidad con una carga mas homogenea entre procesos:

- Se excluyo `Homo.faa` (75,8 millones de aminoacidos) por ser
  desproporcionadamente grande.
- Se excluyeron los proteomas de *Wolbachia* (~0,3 millones de aminoacidos)
  por ser desproporcionadamente pequenos.
- Se eliminaron entradas duplicadas presentes en `filelist_55`.
- Se anadieron especies de tamano medio descargadas de UniProt Reference
  Proteomes (https://www.uniprot.org/proteomes) hasta alcanzar 64 especies,
  numero que coincide con el maximo de procesos evaluado en los benchmarks
  (`np = 64`).

Con estos cambios, el ratio entre la especie mas grande y la mas pequena se
reduce de ≈220x (en `filelist_55`) a ≈8,5x. El detalle de las especies
incluidas, su origen y su tamano en aminoacidos se documenta en la memoria
del TFM (capitulo de datasets y metodologia experimental).

No se uso `filelist_40` en ningun informe de avance.
