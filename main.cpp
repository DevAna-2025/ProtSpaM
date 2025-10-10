#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
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

inline double wall_time()
{
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

int main(int argc, char **argv)
{
    double start = wall_time();

    if (argc < 2)
    {
        printHelp();
        exit(EXIT_FAILURE);
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
    printParameters(weight, dc, threshold, patternNumber, output_filename);

    cout << "\n===== ETAPA 1: Carga de patrones =====\n";

    vector<vector<char>> patterns;
    if (!loadPatterns.empty())
    {
        cout << "Cargando patrones desde archivo: " << loadPatterns << endl;
        patterns = parsePatterns(loadPatterns);
    }
    else
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

    cout << "Patrones generados (" << patterns.size() << "):\n";
    for (size_t i = 0; i < patterns.size(); ++i)
    {
        cout << "  Pattern " << i << ": ";
        for (char c : patterns[i]) cout << c;
        cout << endl;
    }

    if (savePatterns)
    {
        ofstream patOut("patterns.txt");
        for (auto &p : patterns)
        {
            for (auto &c : p)
                patOut << c;
            patOut << endl;
        }
    }

    cout << "\n===== ETAPA 2: Parsing de secuencias =====\n";
    vector<Species> species;
    if (inFiles.empty() && !input_filename.empty())
    {
        cout << "Archivo Multifasta detectado: " << input_filename << endl;
        parser(input_filename, species);
    }
    else if (!inFiles.empty())
    {
        cout << "Carpeta detectada con " << inFiles.size() << " archivos.\n";
        sw_parser(inFiles, species, patterns);
    }

    cout << "Número de especies: " << species.size() << endl;

    cout << "\n===== ETAPA 3: Calculando spaced-words =====\n";
    double start_sw = wall_time();
    for (unsigned int i = 0; i < species.size(); ++i)
    {
        for (const auto &pattern : patterns)
        {
            spacedWords(species[i], pattern);
        }
    }
    double end_sw = wall_time();
    cout << "Tiempo spaced-words: " << (end_sw - start_sw) << " s\n";

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
                tooDistant = true;
        }
    }

    double end_matches = wall_time();
    cout << "Tiempo matches: " << (end_matches - start_matches) << " s\n";

    cout << "\n===== ETAPA 5: Guardando matriz de distancias =====\n";
    outputDistanceMatrix(species, output_filename, distance, tooDistant);

    double end_total = wall_time();
    cout << "\n===== RESUMEN =====\n";
    cout << "Tiempo total: " << (end_total - start) << " s\n";

    return 0;
}
