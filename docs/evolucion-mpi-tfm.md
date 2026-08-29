# De OpenMP a MPI: seis versiones de Prot-SpaM

Documento de apoyo a la defensa del TFM. Recorre rama a rama la paralelización
de memoria distribuida de Prot-SpaM, con el código concreto que cambia en cada
salto y por qué cambia.

Versión web navegable (con comparador interactivo de ramas, plegado de bloques
de código y diagramas): ver el artifact enlazado desde la sesión de trabajo.

---

## 00 · El algoritmo y sus cinco etapas

Prot-SpaM estima distancias filogenéticas entre proteomas completos sin
alinearlos, mediante *spaced-word matches*. Un patrón es una máscara binaria
(`1000010000…1001`): las posiciones `1` son *match positions* y las `0` son
*don't-care*, que se puntúan con BLOSUM62 pero pueden diferir. Con peso 6 y 40
*don't-care*, cada patrón mide 46 posiciones.

```cpp
class Word {
  public:
    unsigned long long key;   // 6 aminoácidos × 5 bits = 30 bits
    unsigned int       pos;   // posición de inicio en seq
    bool operator< (const Word& w) const { return key < w.key; }
};

class Species {
  public:
    std::string                   header;        // nombre (10 chars, PHYLIP)
    std::vector<char>             seq;           // proteoma concatenado, 0..23
    std::vector<std::vector<Word>> sorted_words; // [patrón][palabras ordenadas]
    std::vector<int>              starts;        // offsets de cada proteína
};
```

| # | Etapa | Función | Coste |
|---|-------|---------|-------|
| 1 | Generar o cargar patrones | `rasbhari` / `parsePatterns` | despreciable |
| 2 | Leer FASTA y codificar | `sw_parser` / `parser` | E/S, lineal |
| 3 | **Calcular spaced-words** | `spacedWords` | O(N·L·m), memoria alta |
| 4 | **Calcular matches y distancias** | `calc_matches` | **O(N²·L·m)** — domina |
| 5 | Escribir matriz PHYLIP | `outputDistanceMatrix` | despreciable |

La «Fase 3» y la «Fase 4» que dan nombre a las ramas son las etapas 3 y 4.

---

## 01 · `master` → `feat/seq`: desmontar OpenMP

No se puede medir un *speedup* contra una versión ya paralelizada. Se construye
una referencia estrictamente secuencial y reproducible.

### 1.1 Fuera OpenMP, dentro `std::chrono`

```cpp
- #include <omp.h>
+ #include <chrono>

+ inline double wall_time()
+ {
+     return std::chrono::duration<double>(
+                std::chrono::steady_clock::now().time_since_epoch())
+         .count();
+ }
```

`steady_clock` es monótono: no salta si el sistema ajusta la hora durante una
ejecución de horas.

### 1.2 Desaparecen los tres `#pragma`

```cpp
  for (unsigned int i = 0; i < species.size(); ++i) {
      distance[i][i] = 0;
-     #pragma omp parallel for
      for (auto j = species.size() - 1; j > i; --j) {
          double mismatch_rate = calc_matches(species[i], species[j],
                                              weight, dc, threshold,
                                              patterns, outputScores);
          distance[i][j] = calc_distance(mismatch_rate);
          distance[j][i] = distance[i][j];
          if (mismatch_rate > 0.8541) tooDistant = true;
      }
  }
```

El bucle interior sigue recorriendo `j` hacia atrás. Con OpenMP eso ayudaba al
reparto; sin OpenMP es irrelevante, pero **se conserva para que el orden de
acumulación en punto flotante sea idéntico al del original** y las matrices se
puedan comparar byte a byte.

### 1.3 El parámetro `-t` se retira

```cpp
- void parseParameters(int argc, char *argv[], int& weight, int& dc, int& threshold,
-                      int& patterns, int& threads, std::vector<std::string>& inputFileNames,
-                      std::string& output, bool& sp, std::string& lp, bool& outputScores);
+ void parseParameters(int argc, char *argv[], int& weight, int& dc, int& threshold,
+                      int& patterns, std::vector<std::string>& inputFileNames,
+                      std::string& output, bool& sp, std::string& lp, bool& outputScores);
```

```make
- CC = g++ -fopenmp
+ CC = g++
```

### 1.4 Lectura robusta del *filelist*

```cpp
- while (!infile.eof()) {
-     std::getline(infile, line, '\n');
-     inputFileNames.push_back(line);   // ← mete "" al final
- }
+ while (std::getline(infile, line)) {
+     if (!line.empty() && line.back() == '\r') {   // filelists creados en Windows
+         line.pop_back();
+     }
+     if (!line.empty()) {
+         inputFileNames.push_back(line);
+     }
+ }
```

> **Hallazgo verificado.** Esa cadena vacía extra no era inocua: compensaba un
> fallo en `sw_parser`, que recorre `i < fileNames.size() - 1` y descarta
> siempre la última entrada. Al limpiar el parser sin tocar `sw_parser`,
> `feat/seq` pasa a descartar la última especie real. Ver §08.1.

> **Discrepancia README/código.** El README de `feat/seq` anuncia una rutina
> `write_words_checksum` que no existe en el código de la rama.

---

## 02 · `feat/seq` → `fase 3a`: la etapa 3 en MPI

Se paraleliza únicamente el cálculo de *spaced-words*; las etapas 1, 2, 4 y 5
siguen en rank 0. Paso conservador: aísla una etapa y permite medirla.

### 2.1 El Makefile detecta el *wrapper* MPI

```make
MPICXX ?= $(firstword $(shell command -v mpicxx 2>/dev/null \
                          || command -v mpic++ 2>/dev/null \
                          || command -v mpiCC  2>/dev/null))
ifeq ($(strip $(MPICXX)),)
$(error No MPI C++ compiler wrapper found. Load an MPI module, \
        for example: module load openmpi or module load mpich)
endif
CC = $(MPICXX)
CFLAGS = -c -Wall -std=c++11 -I $(IDIR) -O3
```

### 2.2 Reparto estático por bloques contiguos

Puramente aritmético: sin comunicación, cada proceso calcula su rango con la
misma fórmula y llega a la misma respuesta.

```cpp
struct Range { int begin; int count; };

Range species_range(int total, int rank, int size)
{
    int base      = total / size;
    int remainder = total % size;

    Range range;
    range.count = base + (rank < remainder ? 1 : 0);
    range.begin = rank * base + min(rank, remainder);
    return range;
}
```

> **Punto débil que el TFM mide después.** Reparte *número de especies*, no
> *volumen de trabajo*. En `filelist_55`, `Homo.faa` (75,8 M aa) es ~220× mayor
> que los proteomas de *Wolbachia* (~0,3 M aa). De ahí nacen `filelist_64`
> (ratio ≈8,5×) y `filelist_300_desbalanceado`.

### 2.3 Serialización manual de `Species`

MPI solo mueve bloques contiguos de tipos conocidos. Protocolo de dos mensajes:
primero el tamaño, después los datos.

```cpp
void send_char_vector(const vector<char> &values, int dest, int tag)
{
    int size = static_cast<int>(values.size());
    MPI_Send(&size, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);      // 1) cuántos

    if (size > 0)
        MPI_Send(values.data(), size, MPI_CHAR, dest, tag, MPI_COMM_WORLD);  // 2) qué
}

vector<char> recv_char_vector(int source, int tag)
{
    int size = 0;
    MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<char> values(size);                                   // ya sé cuánto reservar
    if (size > 0)
        MPI_Recv(values.data(), size, MPI_CHAR, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return values;
}
```

Las `Word` se transponen de array de estructuras a estructura de arrays: dos
mensajes grandes en lugar de millones de pequeños.

```cpp
void send_words(const vector<Word> &words, int dest)
{
    int word_count = static_cast<int>(words.size());
    MPI_Send(&word_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);
    if (word_count == 0) return;

    vector<unsigned long long> keys(word_count);      // dos arrays homogéneos…
    vector<unsigned int>       positions(word_count);

    for (int i = 0; i < word_count; ++i) {
        keys[i]      = words[i].key;
        positions[i] = words[i].pos;
    }

    MPI_Send(keys.data(),      word_count, MPI_UNSIGNED_LONG_LONG, dest, TAG_WORDS, MPI_COMM_WORLD);
    MPI_Send(positions.data(), word_count, MPI_UNSIGNED,           dest, TAG_WORDS, MPI_COMM_WORLD);
}                                                     // …dos mensajes, no 2·10⁶
```

Alternativa canónica: `MPI_Type_create_struct`. La transposición manual evita
*padding* y alineación, garantiza portabilidad entre implementaciones de MPI y
produce búferes contiguos. Cuesta una copia extra; a cambio, la red va llena.

### 2.4 Patrones idénticos en todos los procesos

`rasbhari` es probabilístico: si cada rank generase los suyos, el resultado no
tendría sentido. Se generan en rank 0 y se difunden.

```cpp
void broadcast_patterns(vector<vector<char>> &patterns, int rank)
{
    int pattern_count = static_cast<int>(patterns.size());
    MPI_Bcast(&pattern_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) patterns.resize(pattern_count);

    for (int i = 0; i < pattern_count; ++i) {
        int pattern_size = rank == 0 ? static_cast<int>(patterns[i].size()) : 0;
        MPI_Bcast(&pattern_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0) patterns[i].resize(pattern_size);
        if (pattern_size > 0)
            MPI_Bcast(patterns[i].data(), pattern_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
}
```

Por eso los *benchmarks* usan `-p patterns_clean.txt`: con patrones fijos, dos
ejecuciones con distinto número de procesos son comparables.

### 2.5 El cómputo local no comunica nada

```cpp
// Fase 3 paralela: cada proceso calcula los spaced-words de su subconjunto local.
// No hay comunicacion dentro de este bucle.
void calculate_local_spaced_words(vector<Species> &local_species,
                                  const vector<vector<char>> &patterns)
{
    for (Species &species : local_species) {
        species.sorted_words.clear();
        for (const vector<char> &pattern : patterns)
            spacedWords(species, pattern);
    }
}
```

### 2.6 La recolección, y el techo de Amdahl

Cada trabajador devuelve sus `sorted_words` etiquetados con el índice global,
de forma que el orden de la matriz no depende del reparto.

```cpp
void recv_species_words(vector<Species> &species, int source)
{
    int local_count = 0;
    MPI_Recv(&local_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < local_count; ++i) {
        int global_index = 0, pattern_count = 0;
        MPI_Recv(&global_index,  1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&pattern_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        vector<vector<Word>> sorted_words(pattern_count);
        for (int p = 0; p < pattern_count; ++p)
            sorted_words[p] = recv_words(source);

        species[global_index].sorted_words = sorted_words;   // posición global
    }
}
```

```cpp
// La paralelizacion termina aqui: los procesos trabajadores ya enviaron sus
// spaced-words a rank 0 y no participan en matches ni escritura de matriz.
if (rank != 0) {
    MPI_Finalize();
    return EXIT_SUCCESS;
}
```

> **Esta es la tesis de la memoria.** La fase 3a paraleliza una etapa O(N·L·m) y
> deja secuencial la etapa O(N²·L·m). Por Amdahl, el *speedup* está acotado por
> `1/(1−f)`; como la etapa 4 domina y crece con N², *f* **disminuye** al crecer
> el problema. Además la recolección introduce un coste que no existía en la
> versión secuencial. La fase 3a no es el resultado del TFM: es la demostración
> experimental de que hay que atacar la fase 4.

---

## 03 · `fase 3a` → `fase 3b`: quién lee el disco

Misma etapa paralelizada, misma descomposición, mismo resultado. Cambia una
sola decisión: en 3a rank 0 lee todo y reparte por la red; en 3b cada proceso
abre directamente los ficheros que le tocan.

### 3.1 Cada rank recorta su trozo del *filelist*

```cpp
else if (!inFiles.empty())
{
    total_species = static_cast<int>(inFiles.size());
    species_range(total_species, rank, size, local_begin, local_count);

    vector<string> local_files;
    local_files.reserve(local_count);
    for (int i = 0; i < local_count; ++i)
        local_files.push_back(inFiles[local_begin + i]);      // solo lo mío

    if (rank == 0) {
        cout << "Carpeta detectada con " << inFiles.size() << " archivos.\n";
        cout << "Lectura distribuida: cada rank lee su bloque local.\n";
    }

    sw_parser(local_files, local_species, patterns);           // E/S en paralelo
}
```

Desaparecen `distribute_species()` y las llamadas a `send_species`/`recv_species`
en esta ruta. Las secuencias nunca viajan por la red.

> **Aquí está el arreglo que `feat/seq` no incorporó.** Con lectura distribuida,
> el fallo de `sw_parser` dejaría de ser «se pierde la última especie» y pasaría
> a ser «cada rank pierde la última de su bloque»: con 32 procesos, hasta 32
> especies, y el número dependería del número de procesos.

```cpp
- for (unsigned int i = 0; i < fileNames.size() - 1; ++i)
+ for (unsigned int i = 0; i < fileNames.size(); ++i)
```

### 3.2 La recolección tiene que devolver más cosas

En 3a rank 0 ya tenía todas las `Species` completas. En 3b nunca ha visto las
especies ajenas: no conoce sus cabeceras ni sus secuencias.

```cpp
void send_species_results(const vector<Species> &local_species, int global_begin, int dest)
{
    int local_count = static_cast<int>(local_species.size());
    MPI_Send(&local_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

    for (int i = 0; i < local_count; ++i) {
        int global_index = global_begin + i;
        MPI_Send(&global_index, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

        send_species(local_species[i], dest);       // ← header + seq + starts

        int pattern_count = static_cast<int>(local_species[i].sorted_words.size());
        MPI_Send(&pattern_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);
        for (int p = 0; p < pattern_count; ++p) { /* keys + positions */ }
    }
}
```

```cpp
if (!inFiles.empty())
    gather_species_results(species, local_species, total_species, local_begin, rank, size);
else
    gather_spaced_words   (species, local_species, total_species, rank, size);
```

> **Conclusión honesta.** 3b elimina la comunicación de entrada, pero como la
> etapa 4 sigue siendo secuencial en rank 0, la de salida crece. En el balance
> total, 3b mueve el tráfico de sitio más que reducirlo. Su ventaja real es
> otra: la memoria de entrada nunca se concentra en un proceso, y la E/S
> aprovecha el sistema de ficheros paralelo. La lección conjunta de 3a y 3b es
> la misma: mientras la fase 4 no se paralelice, todo camino termina en rank 0.

---

## 04 · `fase 3a` → `metacache`: la etapa 4 en MPI

El salto grande del TFM. Tres problemas a la vez: qué proceso calcula qué par,
cómo consigue los datos que no tiene, y cómo evitar que la memoria explote. La
respuesta a los tres: **invertir el orden de los bucles**.

### 4.1 Streaming patrón a patrón

En vez de completar la etapa 3 antes de empezar la 4, el patrón pasa a ser el
bucle exterior: calcular un patrón → comunicar → acumular en los pares →
descartar → siguiente patrón. El pico de memoria se reduce en un factor *m*.

```cpp
struct StreamMatchState
{
    unsigned int skip             = 0;   // índice de avance en spacedWords2
    int          total_mismatches = 0;
    int          total_dc         = 0;
    long long    mismatch_sum     = 0;
    long long    dc_sum           = 0;
    bool         multi_done       = false;
};
```

```cpp
for (int pattern_index = 0; pattern_index < (int)patterns.size(); ++pattern_index)
{
    vector<int> positions = dontcare_positions(patterns[pattern_index]);
    calculate_local_spaced_words_for_pattern(local_species, patterns[pattern_index]);

    // Estos vectores son temporales. Se calculan para el patron actual y se
    // liberan al empezar el siguiente, que es justo lo que reduce memoria.
    vector<vector<int>> local_run_lengths(local_species.size());
    for (int li = 0; li < (int)local_species.size(); ++li)
        local_run_lengths[li] = key_run_lengths(local_species[li].sorted_words[0]);

    for (int remote_index = 0; remote_index < total_species; ++remote_index)
    {
        /* comunicar words de remote_index, acumular pares */
    }
}
```

```cpp
void calculate_local_spaced_words_for_pattern(vector<Species> &local_species,
                                              const vector<char> &pattern)
{
    // En la version pipeline solo conservo en memoria los spaced-words del
    // patron actual. Al pasar al siguiente patron se descarta el anterior.
    for (Species &species : local_species) {
        species.sorted_words.clear();     // ← el patrón anterior muere aquí
        spacedWords(species, pattern);
    }
}
```

Nótese `sorted_words[0]`: en el pipeline solo hay un patrón vivo.

### 4.2 Quién calcula qué par: propiedad y la regla `i < j`

La matriz es simétrica: solo hace falta el triángulo superior. El par (i, j)
con i < j lo calcula el propietario de *i*.

```cpp
vector<vector<int>> assignments = build_block_assignments(species, size);
owners.assign(total_species, 0);
for (int owner = 0; owner < size; ++owner)
    for (int global_index : assignments[owner])
        owners[global_index] = owner;

MPI_Bcast(owners.data(), total_species, MPI_INT, 0, MPI_COMM_WORLD);
```

### 4.3 La matriz de necesidad remota

```cpp
// Un rank necesita una especie remota j si tiene alguna especie local i con
// i < j. Solo se calculan esos pares para no duplicar trabajo en la matriz.
bool rank_needs_remote_species(int remote_index, const vector<int> &local_indices)
{
    for (int global_i : local_indices)
        if (global_i < remote_index)
            return true;
    return false;
}
```

```cpp
vector<vector<int>> build_remote_need_matrix(int total_species, int size,
                                             const vector<int> &local_indices)
{
    vector<int> local_needs(total_species, 0);
    for (int r = 0; r < total_species; ++r)
        local_needs[r] = rank_needs_remote_species(r, local_indices) ? 1 : 0;

    vector<int> gathered_needs(size * total_species, 0);
    MPI_Allgather(local_needs.data(),    total_species, MPI_INT,
                  gathered_needs.data(), total_species, MPI_INT, MPI_COMM_WORLD);

    // remote_need_matrix[j][rank] == 1  ⟺  ese rank necesita la especie j
    vector<vector<int>> remote_need_matrix(total_species, vector<int>(size, 0));
    for (int r = 0; r < size; ++r)
        for (int j = 0; j < total_species; ++j)
            remote_need_matrix[j][r] = gathered_needs[r * total_species + j];

    return remote_need_matrix;
}
```

Intercambio deliberado: una colectiva barata al principio (~300 KB con 256
procesos y 300 especies) a cambio de eliminar envíos inútiles durante toda la
ejecución. Cada envío evitado es tráfico que no cruza la red de interconexión.

### 4.4 La caché de metadatos — el nombre de la rama

Al comparar dos especies hacen falta dos cosas distintas:

- Los `sorted_words` del patrón actual — **cambian en cada iteración**.
- El `header`, la `seq` y los `starts` — **idénticos en todos los patrones**.

Enviar la secuencia completa cinco veces multiplicaría por cinco el tráfico más
pesado. La metadata se comunica una sola vez, antes del bucle de patrones.

```cpp
// La metadata de una especie no cambia entre patrones. En esta variante la
// comunicacion de header, secuencia y starts se hace una sola vez antes del
// pipeline; dentro del bucle de patrones solo se envian los spaced-words.
vector<Species> build_remote_metadata_cache(const vector<Species> &local_species,
                                            const vector<int> &local_indices,
                                            const vector<int> &owners,
                                            const vector<vector<int>> &remote_need_matrix,
                                            int total_species, int rank, int size)
{
    vector<Species> metadata_cache(total_species);

    for (int remote_index = 0; remote_index < total_species; ++remote_index)
    {
        const vector<int> &rank_needs_remote = remote_need_matrix[remote_index];

        bool any = false;                                    // ¿la quiere alguien?
        for (int need : rank_needs_remote) if (need) { any = true; break; }
        if (!any) continue;                                  // nadie: ni se toca

        int owner = owner_of_species(remote_index, owners);
        Species metadata = stream_species_metadata_for_phase4(
                               local_species, local_indices, remote_index,
                               owner, rank, size, rank_needs_remote);

        if (rank_needs_remote[rank])
            metadata_cache[remote_index] = metadata;
    }

    return metadata_cache;
}
```

Dentro del bucle de patrones solo circula lo que cambia:

```cpp
// En cada iteracion del pipeline se envian unicamente los spaced-words del
// patron actual y solo a los ranks que los necesitan para sus pares.
vector<Word> stream_pattern_words_for_phase4(..., const vector<int> &rank_needs_remote)
{
    if (rank == owner) {
        auto it = find(local_indices.begin(), local_indices.end(), global_index);
        const Species &local = local_species[it - local_indices.begin()];
        const vector<Word> &words = local.sorted_words[0];

        for (int dest = 0; dest < size; ++dest)
            if (dest != owner && rank_needs_remote[dest])
                send_words(words, dest);        // ← MPI_Send: uno detrás de otro

        return words;
    }

    if (!rank_needs_remote[rank]) return vector<Word>();
    return recv_words(owner);
}
```

Esa línea es exactamente lo que la rama `isend` corregirá.

### 4.5 Dos precálculos por patrón

**Primero**, las posiciones *don't-care* son fijas para un patrón:

```cpp
// Para cada patron guardo solo las posiciones don't-care. En el calculo de
// matches estas son las unicas posiciones que se vuelven a comparar.
static vector<int> dontcare_positions(const vector<char> &pattern)
{
    vector<int> positions;
    positions.reserve(pattern.size());
    for (int i = 0; i < (int)pattern.size(); ++i)
        if (pattern[i] == '0') positions.push_back(i);
    return positions;
}

// Version ligera de scorePattern: recibe las posiciones don't-care ya
// calculadas para no recorrer el patron completo en cada comparacion.
static void score_dontcare_positions(const vector<char> &sequence1, const vector<char> &sequence2,
                                     unsigned int pos1, unsigned int pos2,
                                     const vector<int> &positions, int &score, int &mismatches)
{
    score = 0; mismatches = 0;
    for (int pat : positions) {                    // 40 iteraciones, no 46 con rama
        int aa1 = (int)sequence1[pos1 + pat];
        int aa2 = (int)sequence2[pos2 + pat];
        score += blosum62[aa1][aa2];
        if (aa1 != aa2) ++mismatches;
    }
}
```

**Segundo**, `multiMatch()` reescanea hacia delante en cada llamada. Se
sustituye por una tabla:

```cpp
// Precalculo de bloques con la misma key. Mantiene la misma idea de multiMatch,
// pero evita repetir esa busqueda dentro de los bucles principales.
static vector<int> key_run_lengths(const vector<Word> &words)
{
    vector<int> lengths(words.size(), 1);
    size_t i = 0;
    while (i < words.size()) {
        size_t j = i + 1;
        while (j < words.size() && words[j].key == words[i].key) ++j;
        for (size_t k = i; k < j; ++k)
            lengths[k] = (int)(j - k);       // palabras restantes en el bloque desde k
        i = j;
    }
    return lengths;
}
```

```cpp
- int bl1 = multiMatch(spacedWords1, i);     // vuelve a escanear cada vez
+ int bl1 = run_lengths1[i];                 // O(1)
```

> **Fidelidad numérica: decisión importante y contraintuitiva.** En
> `process_streamed_pattern`, `state.skip`, `state.total_mismatches` y
> `state.total_dc` **no se reinician entre patrones**, y al final de cada patrón
> se hace `state.mismatch_sum += state.total_mismatches`. Parece un error de
> acumulación, y no lo es: el `calc_matches` original declara esas variables
> *fuera* del bucle de patrones y hace `mismatchDontCare.emplace_back(...)` con
> los totales **acumulados**. La versión distribuida reproduce esa peculiaridad
> deliberadamente para que la matriz sea idéntica a la del programa original.
> Si algún día se corrigiese el sesgo, habría que hacerlo en las dos rutas a la vez.

### 4.6 Reducción de la matriz y ahorro de memoria en rank 0

Cada par lo calcula exactamente un proceso y el resto aporta cero: la suma
reconstruye la matriz sin conflictos.

```cpp
vector<double> global_distance;
if (rank == 0) global_distance.resize(total_species * total_species, 0.0);

MPI_Reduce(local_distance.data(),
           rank == 0 ? global_distance.data() : nullptr,
           total_species * total_species, MPI_DOUBLE,
           MPI_SUM, 0, MPI_COMM_WORLD);

int local_flag = local_too_distant ? 1 : 0, global_flag = 0;
MPI_Reduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);  // OR lógico
```

```cpp
if (rank == 0) {
    output_headers = collect_headers(species);   // guardo solo los nombres…
}
distribute_species(species, local_species, local_indices, owners, total_species, rank, size);

if (rank == 0) {
    species.clear();
    species.shrink_to_fit();                     // …y devuelvo la RAM al sistema
}
```

> **La ruta `size == 1`.** Con un solo proceso el programa se desvía a la
> implementación clásica:
> ```cpp
> if (size == 1)
>     return calculate_distance_matrix_sequential(local_species, weight, dc, threshold,
>                                                 patterns, outputScores, tooDistant);
> ```
> Consecuencia directa: **`calc_matches()` solo se ejecuta cuando `np = 1`**. Es
> el motivo por el que existe `isend_opt` y por el que su efecto solo se observa
> en el punto de referencia.

> **Limitación conocida.** La opción `-r` (spamogramas) no está soportada en la
> ruta paralela: `process_streamed_pattern` no mantiene el mapa de *scores*. El
> parámetro se acepta y se ignora cuando `np > 1`.

---

## 05 · `metacache` → `isend`: dejar de esperar en fila

Mismo algoritmo, mismo reparto, mismo resultado numérico. Cambia sólo cómo se
emiten los envíos cuando un propietario manda el mismo dato a varios destinos.

### 5.1 El problema: el protocolo rendezvous

MPI usa dos protocolos según el tamaño. Los mensajes grandes —y los
`sorted_words` de un proteoma lo son— usan **rendezvous**: `MPI_Send` no retorna
hasta que el receptor publica su `MPI_Recv`. Con un bucle de `MPI_Send`
bloqueantes, el propietario se para en el primer destino que no está listo y los
demás esperan turno. El coste es la **suma** de las esperas, no el máximo.

### 5.2 La estructura `PendingSend`

`MPI_Isend` no garantiza haber copiado el búfer al retornar. Si el búfer es
local y se destruye, MPI leerá memoria liberada.

```cpp
// ---------------------------------------------------------------------------
// Envio no bloqueante multi-destino (fix de comunicacion de fase 4).
//
// Problema original: cuando un rank es "owner" de una especie, enviaba sus
// datos a cada destino con MPI_Send bloqueante, uno detras de otro. Con
// mensajes grandes MPI usa rendezvous (el Send no vuelve hasta que el receptor
// posta su Recv), asi que el owner quedaba esperando a cada destino en fila.
//
// Solucion: lanzar todos los Isend de una tanda (mismo dato, varios destinos)
// sin esperar entre ellos, acumular los MPI_Request y hacer un unico Waitall al
// final. Asi el progreso de todos los envios ocurre en paralelo. Los buffers se
// guardan en la estructura PendingSend y deben seguir vivos hasta el Waitall,
// porque Isend no garantiza haber copiado el contenido al retornar.
// ---------------------------------------------------------------------------

struct PendingSend
{
    vector<MPI_Request>        requests;
    vector<int>                size_buffer;
    vector<char>               char_buffer;
    vector<int>                int_buffer;
    vector<unsigned long long> key_buffer;
    vector<unsigned int>       pos_buffer;
};
```

Los búferes se rellenan una sola vez y los N destinos comparten el mismo
puntero. Es correcto: MPI solo *lee* del búfer de envío, y varias peticiones
pueden leer la misma región concurrentemente. Además ahorra N copias del vector
de claves. Lo que sería un error es reutilizar ese búfer antes del `Waitall`.

### 5.3 La primitiva multi-destino

```cpp
void isend_words_multi(const vector<Word> &words, const vector<int> &dests,
                       int tag, PendingSend &pending)
{
    int word_count = static_cast<int>(words.size());
    pending.size_buffer.assign(1, word_count);

    pending.key_buffer.resize(word_count);           // transposición AoS → SoA,
    pending.pos_buffer.resize(word_count);           // UNA vez para todos los destinos
    for (int i = 0; i < word_count; ++i) {
        pending.key_buffer[i] = words[i].key;
        pending.pos_buffer[i] = words[i].pos;
    }

    for (int dest : dests) {                         // el bucle ya no bloquea
        MPI_Request req_count;
        MPI_Isend(pending.size_buffer.data(), 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &req_count);
        pending.requests.push_back(req_count);

        if (word_count > 0) {
            MPI_Request req_keys;
            MPI_Isend(pending.key_buffer.data(), word_count, MPI_UNSIGNED_LONG_LONG,
                      dest, tag, MPI_COMM_WORLD, &req_keys);
            pending.requests.push_back(req_keys);

            MPI_Request req_pos;
            MPI_Isend(pending.pos_buffer.data(), word_count, MPI_UNSIGNED,
                      dest, tag, MPI_COMM_WORLD, &req_pos);
            pending.requests.push_back(req_pos);
        }
    }
}

void wait_all_pending(PendingSend &pending)
{
    if (!pending.requests.empty())
        MPI_Waitall((int)pending.requests.size(), pending.requests.data(), MPI_STATUSES_IGNORE);
}
```

### 5.4 Los dos puntos de llamada que cambian

```cpp
+ vector<int> dests;
  for (int dest = 0; dest < size; ++dest) {
      if (dest != owner && rank_needs_remote[dest]) {
-         send_words(words, dest);          // bloquea aquí, destino a destino
+         dests.push_back(dest);            // solo apunto a quién hay que mandarlo
      }
  }

+ if (!dests.empty()) {
+     PendingSend pending;
+     isend_words_multi(words, dests, TAG_WORDS, pending);   // publica TODOS los Isend
+     wait_all_pending(pending);                            // una sola espera
+ }

  return words;
```

```cpp
+ vector<int> dests;
  for (int dest = 0; dest < size; ++dest) {
      if (dest != owner && rank_needs_remote[dest]) {
-         send_species(metadata, dest);
+         dests.push_back(dest);
      }
  }

+ if (!dests.empty()) {
+     PendingSend pending_header, pending_seq, pending_starts;
+     isend_species_multi(metadata, dests, TAG_SPECIES,
+                         pending_header, pending_seq, pending_starts);
+     wait_all_pending(pending_header);
+     wait_all_pending(pending_seq);
+     wait_all_pending(pending_starts);
+ }
```

> **Matiz que el propio README de `isend_opt` reconoce.** El esquema es no
> bloqueante, pero **no solapa comunicación con cómputo**: entre el `Isend` y el
> `Waitall` no se ejecuta trabajo útil. Lo que se gana es real —el coste pasa de
> suma a máximo— pero no es solapamiento. Un solapamiento auténtico exigiría
> publicar los `Isend`, avanzar los pares que ya tienen datos, y llamar a
> `Waitall` al final de la iteración. Es la extensión natural del trabajo.

**Lo que no cambia entre metacache e isend:** verificado con `git diff`, solo se
tocan `main.cpp` (+164 líneas), el README y `src/calc_matches.cpp`. El reparto,
la matriz de necesidad, la caché, el pipeline y `process_streamed_pattern` son
idénticos. Es un experimento limpio.

---

## 06 · `isend` → `isend_opt`: la optimización que no funcionó

### 6.1 Qué cambia exactamente

| Fichero | `isend` vs `isend_opt` |
|---|---|
| `main.cpp` | **byte a byte idénticos** — el diff está vacío |
| `src/calc_matches.cpp` | 182 líneas modificadas |
| `README.md` | 53 líneas |

Y el `calc_matches.cpp` resultante es **exactamente el de `metacache`**:
`git diff metacache isend_opt -- src/calc_matches.cpp` no devuelve nada. La
historia del fichero es:

```
master  ──►  metacache  ──►  isend  ──►  isend_opt
original     + precálculo    revertido    reintroducido
             run_lengths     = master     = metacache
```

La reversión intermedia no es un descuido: es el diseño experimental. Al volver
`calc_matches.cpp` a su forma original en `isend`, cada rama posterior aísla una
única variable respecto a la anterior.

### 6.2 Dónde se aplica

```cpp
for (unsigned int currentPattern = 0; currentPattern < patterns.size(); ++currentPattern) {
    const vector<Word>& spacedWords1 = species1.sorted_words[currentPattern];
    const vector<Word>& spacedWords2 = species2.sorted_words[currentPattern];
    const vector<char>& pattern      = patterns[currentPattern];
+   const vector<int> positions     = dontcare_positions(pattern);
+   const vector<int> run_lengths1  = key_run_lengths(spacedWords1);
+   const vector<int> run_lengths2  = key_run_lengths(spacedWords2);

    for (unsigned int i = 0; i < spacedWords1.size(); ++i) {
-       int bl1 = multiMatch(spacedWords1, i);
+       int bl1 = run_lengths1[i];
```

Más una limpieza de C++ que no altera el resultado:

```cpp
- vector<int> best = {threshold - 1, 0};   // asignación dinámica en bucle caliente
+ int best_score      = threshold - 1;
+ int best_mismatches = 0;

- for (i; i <= limit1; ++i) {              // expresión sin efecto → warning
+ for (;  i <= limit1; ++i) {
```

### 6.3 El resultado

> **Resultado negativo, conservado como tal.** Según el README de la propia
> rama, con `np = 1`: **`isend_opt` registró 3755,3 s frente a los 3141,4 s de
> `isend`** — una degradación en torno al 20 %.
>
> El README lo dice sin adornos: *«La optimización se evaluó con la hipótesis de
> reducir trabajo repetido en la referencia secuencial. Sin embargo, en los
> ensayos recogidos en el TFM no produjo una mejora […]. Se conserva como
> resultado experimental y como referencia para futuros ajustes»*.

Explicaciones plausibles que conviene llevar preparadas:

1. **Coste de materializar `key_run_lengths`.** Construye un `vector<int>` del
   tamaño completo del array de *spaced-words* —millones de entradas, decenas de
   MB— para cada par y cada patrón. En `calc_matches` ese coste se paga N²/2
   veces. En la ruta paralela, `local_run_lengths` se calcula **una vez por
   especie y patrón** y se reutiliza en todos los pares: ahí sí amortiza. La
   misma optimización en dos contextos distintos, rentable sólo en uno.
2. **Presión sobre la caché.** `multiMatch()` escanea datos ya en caché. La
   tabla precalculada añade un array grande que compite con `seq1`, `seq2` y los
   dos vectores de `Word`.
3. **Presión de asignación.** Tres vectores nuevos por patrón y por par:
   N²/2 × m asignaciones dinámicas grandes.
4. **El compilador ya lo hacía.** Con `-O3`, `multiMatch` es pequeña y estática;
   es razonable que GCC la incorpore.

**Qué mediría antes de la defensa:** una ejecución con `perf stat` sobre
`np = 1` en ambas ramas, comparando *cache-misses*, *instructions* e *IPC*. Si
`isend_opt` ejecuta menos instrucciones pero tiene más fallos de caché, las
hipótesis 1+2 quedan confirmadas con datos.

---

## 07 · Las seis ramas de un vistazo

| Rama | Etapa 2 (lectura) | Etapa 3 | Etapa 4 | Comunicación | Aporta |
|---|---|---|---|---|---|
| `master` | rank único | OpenMP | OpenMP | memoria compartida | punto de partida |
| `feat/seq` | secuencial | secuencial | secuencial | ninguna | referencia de corrección y de tiempos |
| `mpi-phase3-a` | centralizada en rank 0 | **MPI** | rank 0 | `Bcast` + `Send`/`Recv` | primera descomposición; evidencia del techo de Amdahl |
| `mpi-phase3-b` | **distribuida** | **MPI** | rank 0 | `Bcast` + recolección | E/S paralela; memoria de entrada repartida |
| `…metacache` | centralizada | **MPI, streaming** | **MPI** | `Allgather` + `Send` + `Reduce` | fase cuadrática distribuida; pico de memoria ÷ m |
| `…isend` | centralizada | **MPI, streaming** | **MPI** | **`Isend` + `Waitall`** | envíos multi-destino en paralelo: suma → máximo |
| `…isend-calcopt` | centralizada | **MPI, streaming** | **MPI** | `Isend` + `Waitall` | precálculo en la ruta `np = 1`; resultado negativo documentado |

### Mecanismos MPI y qué justifica cada uno

| Primitiva | Dónde | Por qué esa y no otra |
|---|---|---|
| `MPI_Bcast` | patrones, `total_species`, `owners` | Un emisor, todos receptores, dato pequeño e idéntico. Árbol: O(log P) en vez de O(P). |
| `MPI_Send`/`MPI_Recv` | reparto inicial de especies | Cada destino recibe un dato *distinto*; ocurre una sola vez, fuera del bucle caliente. |
| `MPI_Allgather` | matriz de necesidad remota | Todos aportan y todos necesitan el resultado completo. Sustituye P broadcasts. |
| `MPI_Isend` + `MPI_Waitall` | mismo dato a varios destinos | El subconjunto de destinos varía, así que no encaja en ninguna colectiva. Sin bloqueo, coste de suma a máximo. |
| `MPI_Reduce` (`SUM`) | matriz de distancias | Cada par lo calcula un solo proceso y el resto aporta cero: la suma reconstruye la matriz sin carreras. |
| `MPI_Reduce` (`MAX`) | flag `tooDistant` | Con valores 0/1, el máximo *es* el OR lógico. |

---

## 08 · Hallazgos verificados en esta revisión

Comprobados ejecutando el código, no deducidos de la lectura. Ninguno afecta a
la corrección de los resultados de fase 4, que son el núcleo del TFM, pero
varios afectan a la rama de referencia y a las de fase 3.

### 8.1 `feat/seq` descarta la última especie del *filelist*

**Severidad: alta** — afecta a la rama que sirve de patrón de corrección.

```
$ wc -l fl5
5 fl5

$ ./bin/Debug/protspam -w 6 -d 40 -m 5 -l fl5 -p patterns.txt -o DMat5
Número de especies: 4          ← deberían ser 5

$ head -1 DMat5
     4                          ← matriz 4×4; sp5 nunca se leyó

# tras aplicar el mismo arreglo que ya llevan las ramas MPI:
#   - for (unsigned int i = 0; i < fileNames.size() - 1; ++i)
#   + for (unsigned int i = 0; i < fileNames.size(); ++i)
$ ./bin/Debug/protspam -w 6 -d 40 -m 5 -l fl5 -p patterns.txt -o DMat5fix
Número de especies: 5          ✓
```

**Implicación.** Cualquier comparación de tiempos entre `feat/seq` y las ramas
MPI sobre el mismo *filelist* enfrenta *n−1* especies contra *n*. Con
`filelist_300`, 44 551 pares frente a 44 850: en torno a un 0,7 % de trabajo de
menos para la referencia, en la dirección de **infraestimar** el *speedup*.
Además las matrices tienen dimensiones distintas, así que la verificación de
corrección no puede hacerse por comparación directa.

### 8.2 Tres ramas no compilan con GCC 13

**Severidad: media.**

```
$ make
src/calc_matches.cpp:59:23: error: 'INT32_MIN' was not declared in this scope
make: *** [Makefile:29: obj/calc_matches.o] Error 1
```

Las versiones recientes de libstdc++ dejaron de incluir `<cstdint>` de forma
transitiva. El repositorio fue añadiendo ese *include* fichero a fichero, y
`calc_matches.cpp` solo lo recibió a partir de la fase 4:

| Rama | `rasbopt.cpp` | `parameters.cpp` | `calc_matches.cpp` | ¿compila con GCC 13? |
|---|---|---|---|---|
| `master` | — | — | — | no |
| `feat/seq` | sí | sí | — | **no** |
| `mpi-phase3-a` | sí | sí | — | **no** |
| `mpi-phase3-b` | sí | sí | — | **no** |
| `…metacache` | sí | sí | sí | sí |
| `…isend` | sí | sí | sí | sí |
| `…isend-calcopt` | sí | sí | sí | sí |

```cpp
  #include "calc_matches.h"
+ #include <cstdint>
```

### 8.3 El README de `feat/seq` describe una función que no existe

**Severidad: media** — es una afirmación sobre el método de validación.

El README anuncia `write_words_checksum` como el mecanismo con el que se
verifica que las versiones paralelas producen el mismo resultado. La función no
está en el código (`git grep` solo la encuentra en el propio README). O se
implementa —volcar `sorted_words` a fichero y comparar con `diff` da un
argumento de corrección mucho más fuerte que comparar sólo la matriz final— o se
retira la frase.

### 8.4 La opción `-r` se ignora en la ruta paralela

**Severidad: baja** — limitación funcional conocida.

`process_streamed_pattern` no lleva el mapa de *scores*, de modo que los
spamogramas sólo se generan con `np = 1`. Con `np > 1` el parámetro se acepta en
silencio. Bastaría con un aviso en rank 0 y una línea en la memoria.

---

## Fuentes

Todo el contenido procede del análisis directo de las seis ramas de
`ana-izaguirre/ProtSpaM`: `git diff` entre ramas consecutivas, lectura del
código y, en §08, ejecución real del binario. Las cifras de `np = 1` para
`isend` e `isend_opt` (3141,4 s y 3755,3 s) proceden del README de la rama
`feat/mpi-phase4-metacache-isend-calcopt`.

Herramienta base: *Prot-SpaM*, Leimeister et al., *GigaScience* 8(3), giy148
(2019) — <https://github.com/jschellh/ProtSpaM>. Proteomas:
<http://projects.gobics.de/data/protspam/paperData.tgz> y UniProt Reference
Proteomes.
