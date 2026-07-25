# Filelists (Fase 4)

Los *filelists* son listas de texto plano con las rutas a los proteomas usados
en cada experimento. Definen qué especies se incluyen y en qué orden, pero no
contienen las secuencias en sí.

## Origen de los proteomas

Los proteomas provienen, salvo que se indique lo contrario, del conjunto de
datos del artículo original de Prot-SpaM:

http://projects.gobics.de/data/protspam/paperData.tgz

Los proteomas no se incluyen en este repositorio (véase la sección
"Datos de entrada" en el README de la rama). Para reproducir los experimentos,
descargar el `.tgz` anterior y extraerlo en la carpeta `data/`. Los conjuntos
de mayor tamaño se completan además con especies de UniProt Reference Proteomes
(https://www.uniprot.org/proteomes), como se detalla más abajo.

## filelist_55

Conjunto de 55 entradas (54 especies únicas). Incluye `Homo.faa`, la especie de
mayor tamaño del conjunto (75,8 millones de aminoácidos), junto con las de menor
tamaño (proteomas de *Wolbachia*, ~0,3 millones de aminoácidos). Esta diferencia
de tamaño entre especies (ratio máximo/mínimo ≈ 220×) provoca un desbalance de
carga entre procesos, analizado en la memoria del TFM.

## filelist_64

Conjunto de dispersión reducida de 64 especies, construido a partir de
`filelist_55` para evaluar la escalabilidad con una carga más homogénea entre
procesos:

- Se excluyó `Homo.faa` (75,8 millones de aminoácidos) por ser
  desproporcionadamente grande.
- Se excluyeron los proteomas de *Wolbachia* (~0,3 millones de aminoácidos) por
  ser desproporcionadamente pequeños.
- Se eliminaron las entradas duplicadas presentes en `filelist_55`.
- Se añadieron especies de tamaño medio de UniProt Reference Proteomes hasta
  alcanzar 64 especies, número que coincide con el máximo de procesos evaluado
  en los benchmarks de nodo único (`np = 64`).

Con estos cambios, el ratio entre la especie más grande y la más pequeña se
reduce de ≈220× (en `filelist_55`) a ≈8,5×.

## filelist_300 / filelist_300_desbalanceado

Conjuntos de 300 especies, los mayores evaluados, usados en las pruebas de
escalabilidad multinodo (hasta `np = 256`):

- `filelist_300`: conjunto de dispersión reducida, con tamaños de proteoma
  relativamente homogéneos entre especies.
- `filelist_300_desbalanceado`: mismo número de especies, pero con una
  distribución de tamaños deliberadamente dispar, para medir el efecto del
  desbalance de carga entre procesos.

El contraste entre ambos permite aislar la penalización por desbalance de la
propia del esquema de paralelización.

## Documentación

El detalle de las especies incluidas en cada conjunto, su origen y su tamaño en
aminoácidos se documenta en la memoria del TFM (capítulo de datasets y
metodología experimental).
