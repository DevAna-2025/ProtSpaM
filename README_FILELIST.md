# Filelists (Fase 3)

Los *filelists* son listas de texto plano con las rutas a los proteomas usados
en cada experimento. Definen qué especies se incluyen y en qué orden, pero no
contienen las secuencias en sí.

## Origen de los proteomas

Los proteomas provienen del conjunto de datos del artículo original de
Prot-SpaM:

http://projects.gobics.de/data/protspam/paperData.tgz

Los proteomas no se incluyen en este repositorio (véase la sección
"Datos de entrada" en el README de la rama). Para reproducir los experimentos,
descargar el `.tgz` anterior y extraerlo en la carpeta `data/`.

## filelist_10 / filelist_20 / filelist_30

Subconjuntos anidados y crecientes de especies (10, 20 y 30), usados en la
evaluación de la Fase 3 para caracterizar la escalabilidad al aumentar
progresivamente el tamaño del problema. Todos incluyen `Homo.faa`.

## Documentación

El detalle de las especies incluidas en cada conjunto, su origen y su tamaño en
aminoácidos se documenta en la memoria del TFM (capítulo de datasets y
metodología experimental).
