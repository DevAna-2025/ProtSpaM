#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <mpi.h>
#include "Species.h"
#include "misc.h"
#include "calc_matches.h"
#include "parser.h"
#include "rasbimp.hpp"
#include "rasbhari.hpp"
#include "patternset.hpp"
#include "parameters.h"
#include "sw_parser.h"
#include <chrono>

using namespace std;

enum MpiTag
{
    TAG_SPECIES = 100,
    TAG_WORDS = 200
};

inline double wall_time()
{
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

void species_range(int total, int rank, int size, int &begin, int &count)
{
    int base = total / size;
    int remainder = total % size;
    count = base + (rank < remainder ? 1 : 0);
    begin = rank * base + min(rank, remainder);
}

void broadcast_patterns(vector<vector<char>> &patterns, int rank)
{
    int pattern_count = static_cast<int>(patterns.size());
    MPI_Bcast(&pattern_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0)
    {
        patterns.resize(pattern_count);
    }

    for (int p = 0; p < pattern_count; ++p)
    {
        int pattern_size = rank == 0 ? static_cast<int>(patterns[p].size()) : 0;
        MPI_Bcast(&pattern_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank != 0)
        {
            patterns[p].resize(pattern_size);
        }
        if (pattern_size > 0)
        {
            MPI_Bcast(patterns[p].data(), pattern_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        }
    }
}

void send_species(const Species &species, int dest)
{
    int header_size = static_cast<int>(species.header.size());
    MPI_Send(&header_size, 1, MPI_INT, dest, TAG_SPECIES, MPI_COMM_WORLD);
    if (header_size > 0)
    {
        MPI_Send(species.header.data(), header_size, MPI_CHAR, dest, TAG_SPECIES, MPI_COMM_WORLD);
    }

    int seq_size = static_cast<int>(species.seq.size());
    MPI_Send(&seq_size, 1, MPI_INT, dest, TAG_SPECIES, MPI_COMM_WORLD);
    if (seq_size > 0)
    {
        MPI_Send(species.seq.data(), seq_size, MPI_CHAR, dest, TAG_SPECIES, MPI_COMM_WORLD);
    }

    int starts_size = static_cast<int>(species.starts.size());
    MPI_Send(&starts_size, 1, MPI_INT, dest, TAG_SPECIES, MPI_COMM_WORLD);
    if (starts_size > 0)
    {
        MPI_Send(species.starts.data(), starts_size, MPI_INT, dest, TAG_SPECIES, MPI_COMM_WORLD);
    }
}

Species recv_species(int source)
{
    Species species;

    int header_size = 0;
    MPI_Recv(&header_size, 1, MPI_INT, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    string header(header_size, '\0');
    if (header_size > 0)
    {
        MPI_Recv(&header[0], header_size, MPI_CHAR, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    species.set_header(header);

    int seq_size = 0;
    MPI_Recv(&seq_size, 1, MPI_INT, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    vector<char> seq(seq_size);
    if (seq_size > 0)
    {
        MPI_Recv(seq.data(), seq_size, MPI_CHAR, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    for (char amino_acid : seq)
    {
        species.set_seq(amino_acid);
    }

    int starts_size = 0;
    MPI_Recv(&starts_size, 1, MPI_INT, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    vector<int> starts(starts_size);
    if (starts_size > 0)
    {
        MPI_Recv(starts.data(), starts_size, MPI_INT, source, TAG_SPECIES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    species.set_starts(starts);

    return species;
}

void distribute_species(vector<Species> &species, vector<Species> &local_species,
                        int total_species, int rank, int size)
{
    int local_begin = 0;
    int local_count = 0;
    species_range(total_species, rank, size, local_begin, local_count);

    if (rank == 0)
    {
        for (int dest = 1; dest < size; ++dest)
        {
            int dest_begin = 0;
            int dest_count = 0;
            species_range(total_species, dest, size, dest_begin, dest_count);
            for (int i = 0; i < dest_count; ++i)
            {
                send_species(species[dest_begin + i], dest);
            }
        }

        local_species.assign(species.begin() + local_begin,
                             species.begin() + local_begin + local_count);
    }
    else
    {
        local_species.reserve(local_count);
        for (int i = 0; i < local_count; ++i)
        {
            local_species.push_back(recv_species(0));
        }
    }
}

void send_species_words(const vector<Species> &local_species, int global_begin, int dest)
{
    int local_count = static_cast<int>(local_species.size());
    MPI_Send(&local_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

    for (int i = 0; i < local_count; ++i)
    {
        int global_index = global_begin + i;
        MPI_Send(&global_index, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

        int pattern_count = static_cast<int>(local_species[i].sorted_words.size());
        MPI_Send(&pattern_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

        for (int p = 0; p < pattern_count; ++p)
        {
            const vector<Word> &words = local_species[i].sorted_words[p];
            int word_count = static_cast<int>(words.size());
            MPI_Send(&word_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

            if (word_count > 0)
            {
                vector<unsigned long long> keys(word_count);
                vector<unsigned int> positions(word_count);
                for (int w = 0; w < word_count; ++w)
                {
                    keys[w] = words[w].key;
                    positions[w] = words[w].pos;
                }
                MPI_Send(keys.data(), word_count, MPI_UNSIGNED_LONG_LONG, dest, TAG_WORDS, MPI_COMM_WORLD);
                MPI_Send(positions.data(), word_count, MPI_UNSIGNED, dest, TAG_WORDS, MPI_COMM_WORLD);
            }
        }
    }
}

void recv_species_words(vector<Species> &species, int source)
{
    int local_count = 0;
    MPI_Recv(&local_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < local_count; ++i)
    {
        int global_index = 0;
        MPI_Recv(&global_index, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int pattern_count = 0;
        MPI_Recv(&pattern_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        species[global_index].sorted_words.clear();
        species[global_index].sorted_words.resize(pattern_count);

        for (int p = 0; p < pattern_count; ++p)
        {
            int word_count = 0;
            MPI_Recv(&word_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<Word> words(word_count);

            if (word_count > 0)
            {
                vector<unsigned long long> keys(word_count);
                vector<unsigned int> positions(word_count);
                MPI_Recv(keys.data(), word_count, MPI_UNSIGNED_LONG_LONG, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(positions.data(), word_count, MPI_UNSIGNED, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int w = 0; w < word_count; ++w)
                {
                    words[w].set_key(keys[w]);
                    words[w].set_pos(positions[w]);
                }
            }
            species[global_index].sorted_words[p] = words;
        }
    }
}

void send_species_results(const vector<Species> &local_species, int global_begin, int dest)
{
    int local_count = static_cast<int>(local_species.size());
    MPI_Send(&local_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

    for (int i = 0; i < local_count; ++i)
    {
        int global_index = global_begin + i;
        MPI_Send(&global_index, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);
        send_species(local_species[i], dest);

        int pattern_count = static_cast<int>(local_species[i].sorted_words.size());
        MPI_Send(&pattern_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

        for (int p = 0; p < pattern_count; ++p)
        {
            const vector<Word> &words = local_species[i].sorted_words[p];
            int word_count = static_cast<int>(words.size());
            MPI_Send(&word_count, 1, MPI_INT, dest, TAG_WORDS, MPI_COMM_WORLD);

            if (word_count > 0)
            {
                vector<unsigned long long> keys(word_count);
                vector<unsigned int> positions(word_count);
                for (int w = 0; w < word_count; ++w)
                {
                    keys[w] = words[w].key;
                    positions[w] = words[w].pos;
                }
                MPI_Send(keys.data(), word_count, MPI_UNSIGNED_LONG_LONG, dest, TAG_WORDS, MPI_COMM_WORLD);
                MPI_Send(positions.data(), word_count, MPI_UNSIGNED, dest, TAG_WORDS, MPI_COMM_WORLD);
            }
        }
    }
}

void recv_species_results(vector<Species> &species, int source)
{
    int local_count = 0;
    MPI_Recv(&local_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < local_count; ++i)
    {
        int global_index = 0;
        MPI_Recv(&global_index, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        Species received = recv_species(source);

        int pattern_count = 0;
        MPI_Recv(&pattern_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        received.sorted_words.clear();
        received.sorted_words.resize(pattern_count);

        for (int p = 0; p < pattern_count; ++p)
        {
            int word_count = 0;
            MPI_Recv(&word_count, 1, MPI_INT, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<Word> words(word_count);

            if (word_count > 0)
            {
                vector<unsigned long long> keys(word_count);
                vector<unsigned int> positions(word_count);
                MPI_Recv(keys.data(), word_count, MPI_UNSIGNED_LONG_LONG, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(positions.data(), word_count, MPI_UNSIGNED, source, TAG_WORDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int w = 0; w < word_count; ++w)
                {
                    words[w].set_key(keys[w]);
                    words[w].set_pos(positions[w]);
                }
            }
            received.sorted_words[p] = words;
        }

        species[global_index] = received;
    }
}

void gather_species_results(vector<Species> &species, const vector<Species> &local_species,
                            int total_species, int local_begin, int rank, int size)
{
    if (rank == 0)
    {
        species.resize(total_species);
        for (int i = 0; i < static_cast<int>(local_species.size()); ++i)
        {
            species[local_begin + i] = local_species[i];
        }

        for (int source = 1; source < size; ++source)
        {
            recv_species_results(species, source);
        }
    }
    else
    {
        send_species_results(local_species, local_begin, 0);
    }
}

void gather_spaced_words(vector<Species> &species, const vector<Species> &local_species,
                         int total_species, int rank, int size)
{
    int local_begin = 0;
    int local_count = 0;
    species_range(total_species, rank, size, local_begin, local_count);

    if (rank == 0)
    {
        for (int i = 0; i < local_count; ++i)
        {
            species[local_begin + i].sorted_words = local_species[i].sorted_words;
        }

        for (int source = 1; source < size; ++source)
        {
            recv_species_words(species, source);
        }
    }
    else
    {
        send_species_words(local_species, local_begin, 0);
    }
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

    parseParameters(argc, argv, weight, dc, threshold, patternNumber, inFiles, output_filename, savePatterns, loadPatterns, outputScores);
    string input_filename(argv[argc - 1]);

    if (rank == 0)
    {
        printParameters(weight, dc, threshold, patternNumber, output_filename);
        cout << "Modo MPI: opcion B - lectura distribuida por rank para entradas -l.\n";
        cout << "\n===== ETAPA 1: Carga de patrones =====\n";
    }

    vector<vector<char>> patterns;
    if (rank == 0 && !loadPatterns.empty())
    {
        cout << "Cargando patrones desde archivo: " << loadPatterns << endl;
        patterns = parsePatterns(loadPatterns);
    }
    else if (rank == 0)
    {
        cout << "Generando patrones aleatorios (rasbhari)...\n";
        rasbhari rasb_set = rasb_implement::hillclimb_oc(patternNumber, weight, dc, dc);
        patternset PatSet = rasb_set.pattern_set();
        for (const auto &i : PatSet)
        {
            string pattern = i.to_string();
            vector<char> tmp(pattern.begin(), pattern.end());
            patterns.push_back(tmp);
        }
    }

    broadcast_patterns(patterns, rank);

    if (rank == 0)
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

    if (rank == 0 && savePatterns)
    {
        ofstream patOut("patterns.txt");
        for (auto &p : patterns)
        {
            for (auto &c : p)
            {
                patOut << c;
            }
            patOut << endl;
        }
    }

    if (rank == 0)
    {
        cout << "\n===== ETAPA 2: Parsing de secuencias =====\n";
    }

    vector<Species> species;
    vector<Species> local_species;
    int total_species = 0;
    int local_begin = 0;
    int local_count = 0;

    if (inFiles.empty() && !input_filename.empty())
    {
        if (rank == 0)
        {
            cout << "Archivo Multifasta detectado: " << input_filename << endl;
            parser(input_filename, species);
        }

        total_species = rank == 0 ? static_cast<int>(species.size()) : 0;
        MPI_Bcast(&total_species, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            cout << "Numero de especies: " << species.size() << endl;
        }

        distribute_species(species, local_species, total_species, rank, size);
    }
    else if (!inFiles.empty())
    {
        total_species = static_cast<int>(inFiles.size());
        species_range(total_species, rank, size, local_begin, local_count);

        vector<string> local_files;
        local_files.reserve(local_count);
        for (int i = 0; i < local_count; ++i)
        {
            local_files.push_back(inFiles[local_begin + i]);
        }

        if (rank == 0)
        {
            cout << "Carpeta detectada con " << inFiles.size() << " archivos.\n";
            cout << "Lectura distribuida: cada rank lee su bloque local.\n";
            cout << "Numero de especies: " << total_species << endl;
        }

        sw_parser(local_files, local_species, patterns);
    }

    if (rank == 0)
    {
        cout << "\n===== ETAPA 3: Calculando spaced-words con MPI =====\n";
    }

    double start_sw = wall_time();
    for (unsigned int i = 0; i < local_species.size(); ++i)
    {
        local_species[i].sorted_words.clear();
        for (const auto &pattern : patterns)
        {
            spacedWords(local_species[i], pattern);
        }
    }

    if (!inFiles.empty())
    {
        gather_species_results(species, local_species, total_species, local_begin, rank, size);
    }
    else
    {
        gather_spaced_words(species, local_species, total_species, rank, size);
    }

    double end_sw = wall_time();
    if (rank == 0)
    {
        cout << "Tiempo spaced-words: " << (end_sw - start_sw) << " s\n";
    }

    if (rank != 0)
    {
        MPI_Finalize();
        return 0;
    }

    cout << "\n===== ETAPA 4: Calculando matches =====\n";
    double start_matches = wall_time();
    vector<vector<double>> distance(species.size(), vector<double>(species.size()));

    for (unsigned int i = 0; i < species.size(); ++i)
    {
        distance[i][i] = 0;
        for (auto j = species.size() - 1; j > i; --j)
        {
            double mismatch_rate = calc_matches(species[i], species[j], weight, dc, threshold, patterns, outputScores);
            distance[i][j] = calc_distance(mismatch_rate);
            distance[j][i] = distance[i][j];
            if (mismatch_rate > 0.8541)
            {
                tooDistant = true;
            }
        }
    }

    double end_matches = wall_time();
    cout << "Tiempo matches: " << (end_matches - start_matches) << " s\n";

    cout << "\n===== ETAPA 5: Guardando matriz de distancias =====\n";
    outputDistanceMatrix(species, output_filename, distance, tooDistant);

    double end_total = wall_time();
    cout << "\n===== RESUMEN =====\n";
    cout << "Tiempo total: " << (end_total - start) << " s\n";

    MPI_Finalize();
    return 0;
}
