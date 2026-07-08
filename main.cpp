#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <string>
#include <vector>

#include "Species.h"
#include "calc_matches.h"
#include "misc.h"
#include "parameters.h"
#include "parser.h"
#include "patternset.hpp"
#include "rasbhari.hpp"
#include "rasbimp.hpp"
#include "sw_parser.h"

using namespace std;

extern int blosum62[24][24];

enum MpiTag
{
    TAG_SPECIES = 100,
    TAG_WORDS = 200
};

struct SpeciesAssignment
{
    int index;
    size_t size;
};

struct StreamMatchState
{
    unsigned int skip = 0;
    int total_mismatches = 0;
    int total_dc = 0;
    long long mismatch_sum = 0;
    long long dc_sum = 0;
    bool multi_done = false;
};

// Para cada patron guardo solo las posiciones don't-care. En el calculo de
// matches estas son las unicas posiciones que se vuelven a comparar.
static vector<int> dontcare_positions(const vector<char> &pattern)
{
    vector<int> positions;
    positions.reserve(pattern.size());

    for (int i = 0; i < static_cast<int>(pattern.size()); ++i)
    {
        if (pattern[i] == '0')
        {
            positions.push_back(i);
        }
    }

    return positions;
}

// Version ligera de scorePattern: recibe las posiciones don't-care ya
// calculadas para no recorrer el patron completo en cada comparacion.
static void score_dontcare_positions(const vector<char> &sequence1,
                                     const vector<char> &sequence2,
                                     unsigned int pos1,
                                     unsigned int pos2,
                                     const vector<int> &positions,
                                     int &score,
                                     int &mismatches)
{
    score = 0;
    mismatches = 0;

    for (int pat : positions)
    {
        int aa1 = static_cast<int>(sequence1[pos1 + pat]);
        int aa2 = static_cast<int>(sequence2[pos2 + pat]);
        score += blosum62[aa1][aa2];

        if (aa1 != aa2)
        {
            ++mismatches;
        }
    }
}

// Los spaced-words estan ordenados por key. Esta tabla indica cuantos elementos
// quedan en el bloque actual con la misma key y evita llamar muchas veces a
// multiMatch durante los bucles mas internos.
static vector<int> key_run_lengths(const vector<Word> &words)
{
    vector<int> lengths(words.size(), 1);
    size_t i = 0;

    while (i < words.size())
    {
        size_t j = i + 1;
        while (j < words.size() && words[j].key == words[i].key)
        {
            ++j;
        }

        for (size_t k = i; k < j; ++k)
        {
            lengths[k] = static_cast<int>(j - k);
        }

        i = j;
    }

    return lengths;
}

double wall_time()
{
    return chrono::duration<double>(
               chrono::steady_clock::now().time_since_epoch())
        .count();
}

void send_string(const string &value, int dest, int tag)
{
    int size = static_cast<int>(value.size());
    MPI_Send(&size, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

    if (size > 0)
    {
        MPI_Send(value.data(), size, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    }
}

string recv_string(int source, int tag)
{
    int size = 0;
    MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    string value(size, '\0');
    if (size > 0)
    {
        MPI_Recv(&value[0], size, MPI_CHAR, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return value;
}

void send_char_vector(const vector<char> &values, int dest, int tag)
{
    int size = static_cast<int>(values.size());
    MPI_Send(&size, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

    if (size > 0)
    {
        MPI_Send(values.data(), size, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    }
}

vector<char> recv_char_vector(int source, int tag)
{
    int size = 0;
    MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<char> values(size);
    if (size > 0)
    {
        MPI_Recv(values.data(), size, MPI_CHAR, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return values;
}

void send_int_vector(const vector<int> &values, int dest, int tag)
{
    int size = static_cast<int>(values.size());
    MPI_Send(&size, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

    if (size > 0)
    {
        MPI_Send(values.data(), size, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
}

vector<int> recv_int_vector(int source, int tag)
{
    int size = 0;
    MPI_Recv(&size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<int> values(size);
    if (size > 0)
    {
        MPI_Recv(values.data(), size, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return values;
}

void send_words(const vector<Word> &words, int dest)
{
    int word_count = static_cast<int>(words.size());
    MPI_Send(&word_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

    if (word_count == 0)
    {
        return;
    }

    vector<unsigned long long> keys(word_count);
    vector<unsigned int> positions(word_count);

    for (int i = 0; i < word_count; ++i)
    {
        keys[i] = words[i].key;
        positions[i] = words[i].pos;
    }

    MPI_Send(keys.data(), word_count, MPI_UNSIGNED_LONG_LONG, dest, TAG_WORDS, MPI_COMM_WORLD);
    MPI_Send(positions.data(), word_count, MPI_UNSIGNED, dest, TAG_WORDS, MPI_COMM_WORLD);
}

vector<Word> recv_words(int source)
{
    int word_count = 0;
    MPI_Recv(&word_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<Word> words(word_count);
    if (word_count == 0)
    {
        return words;
    }

    vector<unsigned long long> keys(word_count);
    vector<unsigned int> positions(word_count);

    MPI_Recv(keys.data(), word_count, MPI_UNSIGNED_LONG_LONG, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(positions.data(), word_count, MPI_UNSIGNED, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < word_count; ++i)
    {
        words[i].set_key(keys[i]);
        words[i].set_pos(positions[i]);
    }

    return words;
}

// Los patrones se generan o cargan en rank 0 y se difunden al resto. Asi todos
// los procesos trabajan con la misma configuracion, independientemente del nodo
// en el que se ejecuten.
void broadcast_patterns(vector<vector<char>> &patterns, int rank)
{
    int pattern_count = static_cast<int>(patterns.size());
    MPI_Bcast(&pattern_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0)
    {
        patterns.resize(pattern_count);
    }

    for (int i = 0; i < pattern_count; ++i)
    {
        int pattern_size = rank == 0 ? static_cast<int>(patterns[i].size()) : 0;
        MPI_Bcast(&pattern_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0)
        {
            patterns[i].resize(pattern_size);
        }

        if (pattern_size > 0)
        {
            MPI_Bcast(patterns[i].data(), pattern_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        }
    }
}

void send_species(const Species &species, int dest)
{
    send_string(species.header, dest, TAG_SPECIES);
    send_char_vector(species.seq, dest, TAG_SPECIES);
    send_int_vector(species.starts, dest, TAG_SPECIES);
}

Species recv_species(int source)
{
    Species species;

    species.set_header(recv_string(source, TAG_SPECIES));

    vector<char> seq = recv_char_vector(source, TAG_SPECIES);
    for (char amino_acid : seq)
    {
        species.set_seq(amino_acid);
    }

    vector<int> starts = recv_int_vector(source, TAG_SPECIES);
    species.set_starts(starts);

    return species;
}

// En la opcion A mantengo el parsing en rank 0. Para que el reparto de trabajo
// no dependa del orden del filelist, asigno primero las especies mas grandes al
// proceso que tenga menos carga acumulada.
vector<vector<int>> build_balanced_assignments(const vector<Species> &species, int size)
{
    vector<SpeciesAssignment> items;
    items.reserve(species.size());

    for (int i = 0; i < static_cast<int>(species.size()); ++i)
    {
        SpeciesAssignment item;
        item.index = i;
        item.size = species[i].seq.size();
        items.push_back(item);
    }

    sort(items.begin(), items.end(),
         [](const SpeciesAssignment &a, const SpeciesAssignment &b) {
             return a.size > b.size;
         });

    vector<vector<int>> assignments(size);
    vector<size_t> loads(size, 0);

    for (const SpeciesAssignment &item : items)
    {
        int target = static_cast<int>(min_element(loads.begin(), loads.end()) - loads.begin());
        assignments[target].push_back(item.index);
        loads[target] += item.size;
    }

    for (vector<int> &assignment : assignments)
    {
        sort(assignment.begin(), assignment.end());
    }

    return assignments;
}

void send_int_list(const vector<int> &values, int dest)
{
    int count = static_cast<int>(values.size());
    MPI_Send(&count, 1, MPI_INT, dest, TAG_SPECIES, MPI_COMM_WORLD);

    if (count > 0)
    {
        MPI_Send(values.data(), count, MPI_INT, dest, TAG_SPECIES, MPI_COMM_WORLD);
    }
}

vector<int> recv_int_list(int source)
{
    int count = 0;
    MPI_Recv(&count, 1, MPI_INT, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<int> values(count);
    if (count > 0)
    {
        MPI_Recv(values.data(), count, MPI_INT, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return values;
}

// Rank 0 conserva la lectura/parsing secuencial de entrada y reparte especies a
// los demas procesos usando una asignacion balanceada por longitud de secuencia.
void distribute_species(vector<Species> &species,
                        vector<Species> &local_species,
                        vector<int> &local_indices,
                        vector<int> &owners,
                        int total_species,
                        int rank,
                        int size)
{
    if (rank == 0)
    {
        vector<vector<int>> assignments = build_balanced_assignments(species, size);
        owners.assign(total_species, 0);
        for (int owner = 0; owner < size; ++owner)
        {
            for (int global_index : assignments[owner])
            {
                owners[global_index] = owner;
            }
        }

        for (int dest = 1; dest < size; ++dest)
        {
            send_int_list(assignments[dest], dest);
            for (int global_index : assignments[dest])
            {
                send_species(species[global_index], dest);
            }
        }

        local_indices = assignments[0];
        local_species.reserve(local_indices.size());
        for (int global_index : local_indices)
        {
            local_species.push_back(species[global_index]);
        }
    }
    else
    {
        local_indices = recv_int_list(0);
        local_species.reserve(local_indices.size());
        for (size_t i = 0; i < local_indices.size(); ++i)
        {
            local_species.push_back(recv_species(0));
        }
    }

    if (rank != 0)
    {
        owners.resize(total_species);
    }
    MPI_Bcast(owners.data(), total_species, MPI_INT, 0, MPI_COMM_WORLD);
}

// Fase 3 paralela tradicional: cada proceso calcula todos los patrones de sus
// especies locales. Se usa solo en ejecuciones con un proceso, donde no necesito
// el pipeline por patron.
void calculate_local_spaced_words(vector<Species> &local_species,
                                  const vector<vector<char>> &patterns)
{
    for (Species &species : local_species)
    {
        species.sorted_words.clear();

        for (const vector<char> &pattern : patterns)
        {
            spacedWords(species, pattern);
        }
    }
}

void calculate_local_spaced_words_for_pattern(vector<Species> &local_species,
                                              const vector<char> &pattern)
{
    // En la version pipeline solo conservo en memoria los spaced-words del
    // patron actual. Al pasar al siguiente patron se descarta el anterior.
    for (Species &species : local_species)
    {
        species.sorted_words.clear();
        spacedWords(species, pattern);
    }
}

// Etapa 1 secuencial: los patrones se deciden en rank 0. Mas abajo se difunden
// con MPI_Bcast para que todos los procesos usen la misma configuracion.
vector<vector<char>> load_or_generate_patterns(const string &loadPatterns,
                                               int patternNumber,
                                               int weight,
                                               int dc,
                                               int rank)
{
    vector<vector<char>> patterns;

    if (rank != 0)
    {
        return patterns;
    }

    if (!loadPatterns.empty())
    {
        cout << "Cargando patrones desde archivo: " << loadPatterns << endl;
        return parsePatterns(loadPatterns);
    }

    cout << "Generando patrones aleatorios (rasbhari)...\n";
    rasbhari rasb_set = rasb_implement::hillclimb_oc(patternNumber, weight, dc, dc);
    patternset PatSet = rasb_set.pattern_set();

    for (const auto &pattern : PatSet)
    {
        string pattern_string = pattern.to_string();
        patterns.push_back(vector<char>(pattern_string.begin(), pattern_string.end()));
    }

    return patterns;
}

void print_patterns_used(const vector<vector<char>> &patterns)
{
    cout << "Patrones generados (" << patterns.size() << "):\n";

    for (size_t i = 0; i < patterns.size(); ++i)
    {
        cout << "  Pattern " << i << ": ";
        for (char c : patterns[i])
        {
            cout << c;
        }
        cout << endl;
    }
}

void save_patterns(const vector<vector<char>> &patterns)
{
    ofstream patOut("patterns.txt");

    for (const vector<char> &pattern : patterns)
    {
        for (char c : pattern)
        {
            patOut << c;
        }
        patOut << endl;
    }
}

// Etapa 2 secuencial: en esta primera version MPI, rank 0 mantiene la lectura de
// las secuencias y luego distribuye el trabajo de la etapa 3.
void load_species(vector<Species> &species,
                  vector<string> &inFiles,
                  const string &input_filename,
                  const vector<vector<char>> &patterns)
{
    if (inFiles.empty() && !input_filename.empty())
    {
        cout << "Archivo Multifasta detectado: " << input_filename << endl;
        parser(input_filename, species);
        return;
    }

    if (!inFiles.empty())
    {
        cout << "Carpeta detectada con " << inFiles.size() << " archivos.\n";
        sw_parser(inFiles, species, patterns);
    }
}

vector<string> collect_headers(const vector<Species> &species)
{
    vector<string> headers;
    headers.reserve(species.size());

    for (const Species &item : species)
    {
        headers.push_back(item.header);
    }

    return headers;
}

vector<Species> build_output_species(const vector<string> &headers)
{
    vector<Species> output_species(headers.size());

    for (size_t i = 0; i < headers.size(); ++i)
    {
        output_species[i].set_header(headers[i]);
    }

    return output_species;
}

Species clone_species_metadata(const Species &source)
{
    Species metadata;
    metadata.set_header(source.header);

    for (char amino_acid : source.seq)
    {
        metadata.set_seq(amino_acid);
    }

    vector<int> starts = source.starts;
    metadata.set_starts(starts);

    return metadata;
}

int owner_of_species(int global_index, const vector<int> &owners)
{
    return owners[global_index];
}

// Un rank necesita una especie remota j si tiene alguna especie local i con
// i < j. Solo se calculan esos pares para no duplicar trabajo en la matriz.
bool rank_needs_remote_species(int remote_index, const vector<int> &local_indices)
{
    for (int global_i : local_indices)
    {
        if (global_i < remote_index)
        {
            return true;
        }
    }

    return false;
}

// Matriz pequena de control: remote_need_matrix[j][rank] indica si ese rank
// necesita recibir la especie j durante la fase 4. Esto reduce envios
// innecesarios, sobre todo cuando los procesos estan repartidos entre nodos.
vector<vector<int>> build_remote_need_matrix(int total_species,
                                             int size,
                                             const vector<int> &local_indices)
{
    vector<int> local_needs(total_species, 0);
    for (int remote_index = 0; remote_index < total_species; ++remote_index)
    {
        local_needs[remote_index] = rank_needs_remote_species(remote_index, local_indices) ? 1 : 0;
    }

    vector<int> gathered_needs(size * total_species, 0);
    MPI_Allgather(local_needs.data(),
                  total_species,
                  MPI_INT,
                  gathered_needs.data(),
                  total_species,
                  MPI_INT,
                  MPI_COMM_WORLD);

    vector<vector<int>> remote_need_matrix(total_species, vector<int>(size, 0));
    for (int rank_id = 0; rank_id < size; ++rank_id)
    {
        for (int remote_index = 0; remote_index < total_species; ++remote_index)
        {
            remote_need_matrix[remote_index][rank_id] =
                gathered_needs[rank_id * total_species + remote_index];
        }
    }

    return remote_need_matrix;
}

vector<vector<double>> calculate_distance_matrix_sequential(const vector<Species> &species,
                                                            int weight,
                                                            int dc,
                                                            int threshold,
                                                            const vector<vector<char>> &patterns,
                                                            bool outputScores,
                                                            bool &tooDistant)
{
    vector<vector<double>> distance(species.size(), vector<double>(species.size()));

    for (unsigned int i = 0; i < species.size(); ++i)
    {
        distance[i][i] = 0;

        for (unsigned int j = static_cast<unsigned int>(species.size() - 1); j > i; --j)
        {
            double mismatch_rate = calc_matches(species[i], species[j],
                                                weight, dc, threshold,
                                                patterns, outputScores);
            double dist = calc_distance(mismatch_rate);

            distance[i][j] = dist;
            distance[j][i] = dist;

            if (mismatch_rate > 0.8541)
            {
                tooDistant = true;
            }
        }
    }

    return distance;
}

vector<vector<double>> expand_distance_matrix(const vector<double> &flat_distance,
                                              int total_species)
{
    vector<vector<double>> distance(total_species, vector<double>(total_species));

    for (int i = 0; i < total_species; ++i)
    {
        for (int j = 0; j < total_species; ++j)
        {
            distance[i][j] = flat_distance[i * total_species + j];
        }
    }

    return distance;
}

// Calcula la contribucion de un solo patron para un par de especies. La logica
// sigue de cerca calc_matches original, pero recibe los spaced-words del patron
// actual y acumula el estado del par entre patrones. No reinicio estos contadores
// porque la version secuencial tambien acumula por patron antes de calcular la
// tasa final de mismatches.
void process_streamed_pattern(const Species &species1,
                              const Species &species2,
                              const vector<Word> &spacedWords1,
                              const vector<Word> &spacedWords2,
                              const vector<int> &run_lengths1,
                              const vector<int> &run_lengths2,
                              const vector<int> &positions,
                              int dc,
                              int threshold,
                              StreamMatchState &state)
{
    const vector<char> &sequence1 = species1.seq;
    const vector<char> &sequence2 = species2.seq;

    for (unsigned int i = 0; i < spacedWords1.size(); ++i)
    {
        int bl1 = run_lengths1[i];
        if (bl1 > 1)
        {
            int best_score = threshold - 1;
            int best_mismatches = 0;
            unsigned int limit1 = i + bl1 - 1;
            for (; i <= limit1; ++i)
            {
                bool singleMatch = true;
                for (unsigned int j = state.skip; j < spacedWords2.size() && singleMatch; ++j)
                {
                    int bl2 = run_lengths2[j];
                    if (bl2 > 1)
                    {
                        singleMatch = false;
                        unsigned int limit2 = j + bl2 - 1;
                        for (; j <= limit2; ++j)
                        {
                            int score = 0;
                            int mismatches = 0;
                            if (spacedWords1[i].key > spacedWords2[j].key)
                            {
                                state.skip += bl2;
                                continue;
                            }
                            if (spacedWords1[i].key < spacedWords2[j].key)
                            {
                                break;
                            }
                            if (spacedWords1[i].key == spacedWords2[j].key)
                            {
                                score_dontcare_positions(sequence1, sequence2,
                                                         spacedWords1[i].pos,
                                                         spacedWords2[j].pos,
                                                         positions, score, mismatches);
                            }
                            if (score >= threshold && score > best_score)
                            {
                                best_score = score;
                                best_mismatches = mismatches;
                            }
                            if (i == limit1 && j == limit2)
                            {
                                state.skip += bl2;
                                state.multi_done = true;
                                break;
                            }
                        }
                    }
                    else
                    {
                        if (spacedWords1[i].key > spacedWords2[j].key)
                        {
                            state.skip = j + 1;
                            continue;
                        }
                        if (spacedWords1[i].key < spacedWords2[j].key)
                        {
                            break;
                        }
                        if (spacedWords1[i].key == spacedWords2[j].key)
                        {
                            int score = 0;
                            int mismatches = 0;
                            score_dontcare_positions(sequence1, sequence2,
                                                     spacedWords1[i].pos,
                                                     spacedWords2[j].pos,
                                                     positions, score, mismatches);
                            if (score >= threshold && score > best_score)
                            {
                                best_score = score;
                                best_mismatches = mismatches;
                            }
                            if (i == limit1)
                            {
                                state.skip = j + 1;
                                state.multi_done = true;
                                break;
                            }
                            break;
                        }
                    }
                }
            }
            if (best_score >= threshold)
            {
                state.total_mismatches += best_mismatches;
                state.total_dc += dc;
            }
            if (bl1 > 1 && state.multi_done)
            {
                --i;
                state.multi_done = false;
            }
        }
        else
        {
            bool go_on = true;
            for (unsigned int j = state.skip; j < spacedWords2.size() && go_on; ++j)
            {
                int bl2 = run_lengths2[j];
                if (bl2 > 1)
                {
                    go_on = false;
                    int best_score = threshold - 1;
                    int best_mismatches = 0;
                    unsigned int limit = j + bl2 - 1;
                    for (; j <= limit; ++j)
                    {
                        int score = 0;
                        int mismatches = 0;
                        if (spacedWords1[i].key > spacedWords2[j].key)
                        {
                            state.skip += bl2;
                            continue;
                        }
                        if (spacedWords1[i].key < spacedWords2[j].key)
                        {
                            break;
                        }
                        if (spacedWords1[i].key == spacedWords2[j].key)
                        {
                            score_dontcare_positions(sequence1, sequence2,
                                                     spacedWords1[i].pos,
                                                     spacedWords2[j].pos,
                                                     positions, score, mismatches);
                            if (score >= threshold && score > best_score)
                            {
                                best_score = score;
                                best_mismatches = mismatches;
                            }
                            if (j == limit)
                            {
                                if (best_score > threshold)
                                {
                                    state.total_mismatches += best_mismatches;
                                    state.total_dc += dc;
                                }
                                state.skip += bl2;
                            }
                        }
                    }
                }
                else
                {
                    if (spacedWords1[i].key > spacedWords2[j].key)
                    {
                        state.skip = j + 1;
                        continue;
                    }
                    if (spacedWords1[i].key < spacedWords2[j].key)
                    {
                        break;
                    }
                    if (spacedWords1[i].key == spacedWords2[j].key)
                    {
                        state.skip = j + 1;
                        int score = 0;
                        int mismatches = 0;
                        score_dontcare_positions(sequence1, sequence2,
                                                 spacedWords1[i].pos,
                                                 spacedWords2[j].pos,
                                                 positions, score, mismatches);
                        if (score >= threshold)
                        {
                            state.total_mismatches += mismatches;
                            state.total_dc += dc;
                        }
                        break;
                    }
                }
            }
        }
    }

    state.mismatch_sum += state.total_mismatches;
    state.dc_sum += state.total_dc;
}

// Envia solo metadatos de una especie remota: header, secuencia y starts. Los
// spaced-words se mandan aparte porque cambian en cada patron del pipeline.
Species stream_species_metadata_for_phase4(const vector<Species> &local_species,
                                           const vector<int> &local_indices,
                                           int global_index,
                                           int owner,
                                           int rank,
                                           int size,
                                           const vector<int> &rank_needs_remote)
{
    if (rank == owner)
    {
        auto it = find(local_indices.begin(), local_indices.end(), global_index);
        const Species &local = local_species[it - local_indices.begin()];
        Species metadata = clone_species_metadata(local);
        for (int dest = 0; dest < size; ++dest)
        {
            if (dest != owner && rank_needs_remote[dest])
            {
                send_species(metadata, dest);
            }
        }
        return metadata;
    }

    if (!rank_needs_remote[rank])
    {
        return Species();
    }

    return recv_species(owner);
}

// En cada iteracion del pipeline se envian unicamente los spaced-words del
// patron actual y solo a los ranks que los necesitan para sus pares.
vector<Word> stream_pattern_words_for_phase4(const vector<Species> &local_species,
                                             const vector<int> &local_indices,
                                             int global_index,
                                             int owner,
                                             int rank,
                                             int size,
                                             const vector<int> &rank_needs_remote)
{
    if (rank == owner)
    {
        auto it = find(local_indices.begin(), local_indices.end(), global_index);
        const Species &local = local_species[it - local_indices.begin()];
        const vector<Word> &words = local.sorted_words[0];

        for (int dest = 0; dest < size; ++dest)
        {
            if (dest != owner && rank_needs_remote[dest])
            {
                send_words(words, dest);
            }
        }

        return words;
    }

    if (!rank_needs_remote[rank])
    {
        return vector<Word>();
    }

    return recv_words(owner);
}

// Fase 4 paralela con streaming por patron. La idea principal es no reunir ni
// replicar todos los spaced-words a la vez: se calcula un patron, se comunican
// los datos necesarios, se actualizan los pares de la matriz y se pasa al
// siguiente patron.
vector<vector<double>> calculate_distance_matrix_parallel(vector<Species> &local_species,
                                                          const vector<int> &local_indices,
                                                          const vector<int> &owners,
                                                          int total_species,
                                                          int rank,
                                                          int size,
                                                          int weight,
                                                          int dc,
                                                          int threshold,
                                                          const vector<vector<char>> &patterns,
                                                          bool outputScores,
                                                          bool &tooDistant,
                                                          double &spaced_words_time)
{
    if (size == 1)
    {
        return calculate_distance_matrix_sequential(local_species, weight, dc, threshold,
                                                    patterns, outputScores, tooDistant);
    }

    vector<double> local_distance(total_species * total_species, 0.0);
    bool local_too_distant = false;
    vector<StreamMatchState> pair_states(local_species.size() * total_species);
    vector<vector<int>> remote_need_matrix = build_remote_need_matrix(total_species,
                                                                      size,
                                                                      local_indices);

    for (int pattern_index = 0; pattern_index < static_cast<int>(patterns.size()); ++pattern_index)
    {
        vector<int> positions = dontcare_positions(patterns[pattern_index]);
        double pattern_sw_start = wall_time();
        calculate_local_spaced_words_for_pattern(local_species, patterns[pattern_index]);
        spaced_words_time += wall_time() - pattern_sw_start;

        // Estos vectores son temporales. Se calculan para el patron actual y se
        // liberan al empezar el siguiente, que es justo lo que reduce memoria.
        vector<vector<int>> local_run_lengths(local_species.size());
        for (int local_i = 0; local_i < static_cast<int>(local_species.size()); ++local_i)
        {
            local_run_lengths[local_i] = key_run_lengths(local_species[local_i].sorted_words[0]);
        }

        for (int remote_index = 0; remote_index < total_species; ++remote_index)
        {
            const vector<int> &rank_needs_remote = remote_need_matrix[remote_index];
            bool any_rank_needs_remote = false;
            for (int need : rank_needs_remote)
            {
                if (need)
                {
                    any_rank_needs_remote = true;
                    break;
                }
            }
            if (!any_rank_needs_remote)
            {
                continue;
            }

            // El owner conserva su especie local; los demas solo reciben datos
            // si tienen algun par pendiente con esta especie remota.
            int owner = owner_of_species(remote_index, owners);
            Species remote_species = stream_species_metadata_for_phase4(local_species,
                                                                        local_indices,
                                                                        remote_index,
                                                                        owner,
                                                                        rank,
                                                                        size,
                                                                        rank_needs_remote);
            vector<Word> remote_words = stream_pattern_words_for_phase4(local_species,
                                                                        local_indices,
                                                                        remote_index,
                                                                        owner,
                                                                        rank,
                                                                        size,
                                                                        rank_needs_remote);
            if (!rank_needs_remote[rank])
            {
                continue;
            }

            vector<int> remote_run_lengths = key_run_lengths(remote_words);

            for (int local_i = 0; local_i < static_cast<int>(local_species.size()); ++local_i)
            {
                int global_i = local_indices[local_i];
                if (global_i >= remote_index)
                {
                    continue;
                }

                StreamMatchState &state = pair_states[local_i * total_species + remote_index];
                process_streamed_pattern(local_species[local_i],
                                         remote_species,
                                         local_species[local_i].sorted_words[0],
                                         remote_words,
                                         local_run_lengths[local_i],
                                         remote_run_lengths,
                                         positions,
                                         dc,
                                         threshold,
                                         state);
            }
        }
    }

    for (int remote_index = 0; remote_index < total_species; ++remote_index)
    {
        for (int local_i = 0; local_i < static_cast<int>(local_species.size()); ++local_i)
        {
            int global_i = local_indices[local_i];
            if (global_i >= remote_index)
            {
                continue;
            }

            const StreamMatchState &state = pair_states[local_i * total_species + remote_index];
            double mismatch_rate = state.dc_sum == 0
                                       ? 0.0
                                       : static_cast<double>(state.mismatch_sum) /
                                             static_cast<double>(state.dc_sum);
            double dist = calc_distance(mismatch_rate);

            local_distance[global_i * total_species + remote_index] = dist;
            local_distance[remote_index * total_species + global_i] = dist;

            if (mismatch_rate > 0.8541)
            {
                local_too_distant = true;
            }
        }
    }

    vector<double> global_distance;
    if (rank == 0)
    {
        global_distance.resize(total_species * total_species, 0.0);
    }

    MPI_Reduce(local_distance.data(),
               rank == 0 ? global_distance.data() : nullptr,
               total_species * total_species,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);

    int local_flag = local_too_distant ? 1 : 0;
    int global_flag = 0;
    MPI_Reduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        tooDistant = global_flag != 0;
        return expand_distance_matrix(global_distance, total_species);
    }

    return vector<vector<double>>();
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start = wall_time();

    if (argc < 2)
    {
        if (rank == 0)
        {
            printHelp();
        }

        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int weight = 6;
    int dc = 40;
    int threshold = 0;
    int patternNumber = 5;
    string loadPatterns;
    bool savePatterns = false;
    bool outputScores = false;
    vector<string> inFiles;
    string output_filename = "DMat";
    bool tooDistant = false;

    parseParameters(argc, argv, weight, dc, threshold, patternNumber,
                    inFiles, output_filename, savePatterns, loadPatterns, outputScores);

    if (rank == 0)
    {
        printParameters(weight, dc, threshold, patternNumber, output_filename);
        cout << "\n===== ETAPA 1: Carga de patrones =====\n";
    }

    vector<vector<char>> patterns = load_or_generate_patterns(loadPatterns, patternNumber, weight, dc, rank);
    broadcast_patterns(patterns, rank);

    if (rank == 0)
    {
        print_patterns_used(patterns);
        if (savePatterns)
        {
            save_patterns(patterns);
        }
    }

    if (rank == 0)
    {
        cout << "\n===== ETAPA 2: Parsing de secuencias =====\n";
    }

    vector<Species> species;
    vector<string> output_headers;
    if (rank == 0)
    {
        string input_filename(argv[argc - 1]);
        load_species(species, inFiles, input_filename, patterns);
        output_headers = collect_headers(species);
        cout << "Numero de especies: " << species.size() << endl;
    }

    int total_species = rank == 0 ? static_cast<int>(species.size()) : 0;
    MPI_Bcast(&total_species, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Species> local_species;
    vector<int> local_indices;
    vector<int> owners;
    distribute_species(species, local_species, local_indices, owners, total_species, rank, size);

    if (rank == 0)
    {
        species.clear();
        species.shrink_to_fit();
    }

    if (rank == 0)
    {
        cout << "\n===== ETAPA 3: Calculando spaced-words con MPI =====\n";
    }

    double start_sw = wall_time();
    double spaced_words_time = 0.0;
    if (size == 1)
    {
        calculate_local_spaced_words(local_species, patterns);
        spaced_words_time = wall_time() - start_sw;
    }

    double start_matches = wall_time();
    vector<vector<double>> distance = calculate_distance_matrix_parallel(local_species,
                                                                         local_indices,
                                                                         owners,
                                                                         total_species,
                                                                         rank,
                                                                         size,
                                                                         weight,
                                                                         dc,
                                                                         threshold,
                                                                         patterns,
                                                                         outputScores,
                                                                         tooDistant,
                                                                         spaced_words_time);

    if (rank == 0)
    {
        double phase4_elapsed = wall_time() - start_matches;
        double matches_time = size > 1 ? phase4_elapsed - spaced_words_time : phase4_elapsed;

        cout << "Tiempo spaced-words: " << spaced_words_time << " s\n";
        cout << "\n===== ETAPA 4: Calculando matches con MPI =====\n";
        cout << "Tiempo matches: " << matches_time << " s\n";
    }

    if (rank != 0)
    {
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    cout << "\n===== ETAPA 5: Guardando matriz de distancias =====\n";
    vector<Species> output_species = build_output_species(output_headers);
    outputDistanceMatrix(output_species, output_filename, distance, tooDistant);

    cout << "\n===== RESUMEN =====\n";
    cout << "Tiempo total: " << (wall_time() - start) << " s\n";

    MPI_Finalize();
    return EXIT_SUCCESS;
}
