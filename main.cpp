#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <mpi.h>
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

enum MpiTag
{
    TAG_SPECIES = 100,
    TAG_WORDS = 200
};

struct Range
{
    int begin;
    int count;
};

double wall_time()
{
    return chrono::duration<double>(
               chrono::steady_clock::now().time_since_epoch())
        .count();
}

Range species_range(int total, int rank, int size)
{
    int base = total / size;
    int remainder = total % size;

    Range range;
    range.count = base + (rank < remainder ? 1 : 0);
    range.begin = rank * base + min(rank, remainder);
    return range;
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

// Los patrones se generan o cargan solo en rank 0 y se difunden al resto de
// procesos. De esta forma todos calculan spaced-words con exactamente el mismo
// conjunto de patrones.
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

// Rank 0 conserva la lectura/parsing secuencial de entrada y reparte especies a
// los demas procesos. Este reparto prepara la primera fase paralelizada: la
// etapa 3, calculo de spaced-words.
void distribute_species(vector<Species> &species,
                        vector<Species> &local_species,
                        int total_species,
                        int rank,
                        int size)
{
    Range local = species_range(total_species, rank, size);

    if (rank == 0)
    {
        for (int dest = 1; dest < size; ++dest)
        {
            Range remote = species_range(total_species, dest, size);
            for (int i = 0; i < remote.count; ++i)
            {
                send_species(species[remote.begin + i], dest);
            }
        }

        local_species.assign(species.begin() + local.begin,
                             species.begin() + local.begin + local.count);
        return;
    }

    local_species.reserve(local.count);
    for (int i = 0; i < local.count; ++i)
    {
        local_species.push_back(recv_species(0));
    }
}

// Fase 3 paralela: cada proceso calcula los spaced-words de su subconjunto local.
// No hay comunicacion dentro de este bucle; cada rank trabaja de forma
// independiente sobre las especies que recibio.
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

// Envia a rank 0 los sorted_words ya calculados por un proceso trabajador. Solo
// se transmiten key y pos de cada Word, que son los datos necesarios en fase 4.
void send_species_words(const vector<Species> &local_species, int global_begin, int dest)
{
    int local_count = static_cast<int>(local_species.size());
    MPI_Send(&local_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

    for (int i = 0; i < local_count; ++i)
    {
        int global_index = global_begin + i;
        MPI_Send(&global_index, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

        const vector<vector<Word>> &sorted_words = local_species[i].sorted_words;
        int pattern_count = static_cast<int>(sorted_words.size());
        MPI_Send(&pattern_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

        for (const vector<Word> &words : sorted_words)
        {
            send_words(words, dest);
        }
    }
}

// Rank 0 reconstruye los sorted_words recibidos y los coloca en la posicion
// global de la especie correspondiente para restaurar el orden original.
void recv_species_words(vector<Species> &species, int source)
{
    int local_count = 0;
    MPI_Recv(&local_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < local_count; ++i)
    {
        int global_index = 0;
        int pattern_count = 0;

        MPI_Recv(&global_index, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&pattern_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        vector<vector<Word>> sorted_words(pattern_count);
        for (int p = 0; p < pattern_count; ++p)
        {
            sorted_words[p] = recv_words(source);
        }

        species[global_index].sorted_words = sorted_words;
    }
}

// Cierre de la fase 3 paralela. Rank 0 copia tambien sus propios resultados
// locales porque los calculo sobre local_species, que es una copia de species.
// Despues recibe los resultados del resto de ranks.
void gather_spaced_words(vector<Species> &species,
                         const vector<Species> &local_species,
                         int total_species,
                         int rank,
                         int size)
{
    Range local = species_range(total_species, rank, size);

    if (rank == 0)
    {
        for (int i = 0; i < local.count; ++i)
        {
            species[local.begin + i].sorted_words = local_species[i].sorted_words;
        }

        for (int source = 1; source < size; ++source)
        {
            recv_species_words(species, source);
        }
        return;
    }

    send_species_words(local_species, local.begin, 0);
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

// Etapa 4 secuencial: una vez reunidos todos los sorted_words, solo rank 0
// calcula los matches y la matriz de distancias.
vector<vector<double>> calculate_distance_matrix(const vector<Species> &species,
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
    if (rank == 0)
    {
        string input_filename(argv[argc - 1]);
        load_species(species, inFiles, input_filename, patterns);
        cout << "Numero de especies: " << species.size() << endl;
    }

    int total_species = rank == 0 ? static_cast<int>(species.size()) : 0;
    MPI_Bcast(&total_species, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Species> local_species;
    distribute_species(species, local_species, total_species, rank, size);

    if (rank == 0)
    {
        cout << "\n===== ETAPA 3: Calculando spaced-words con MPI =====\n";
    }

    double start_sw = wall_time();
    calculate_local_spaced_words(local_species, patterns);
    gather_spaced_words(species, local_species, total_species, rank, size);

    if (rank == 0)
    {
        cout << "Tiempo spaced-words: " << (wall_time() - start_sw) << " s\n";
    }

    // La paralelizacion termina aqui: los procesos trabajadores ya enviaron sus
    // spaced-words a rank 0 y no participan en matches ni escritura de matriz.
    if (rank != 0)
    {
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    cout << "\n===== ETAPA 4: Calculando matches =====\n";
    double start_matches = wall_time();
    vector<vector<double>> distance = calculate_distance_matrix(species, weight, dc, threshold,
                                                                patterns, outputScores, tooDistant);
    cout << "Tiempo matches: " << (wall_time() - start_matches) << " s\n";

    cout << "\n===== ETAPA 5: Guardando matriz de distancias =====\n";
    outputDistanceMatrix(species, output_filename, distance, tooDistant);

    cout << "\n===== RESUMEN =====\n";
    cout << "Tiempo total: " << (wall_time() - start) << " s\n";

    MPI_Finalize();
    return EXIT_SUCCESS;
}
