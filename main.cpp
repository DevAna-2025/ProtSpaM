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

    /* -------------------- */
    /* pattern sets */
    vector<vector<char>> patterns;
    if (!loadPatterns.empty())
    {
        patterns = parsePatterns(loadPatterns);
    }
    else
    {
        rasbhari rasb_set = rasb_implement::hillclimb_oc(patternNumber, weight, dc, dc);
        patternset PatSet = rasb_set.pattern_set();
        for (const auto &i : PatSet)
        {
            string pattern = i.to_string();
            vector<char> tmp(pattern.begin(), pattern.end());
            patterns.push_back(tmp);
        }
    }
    print_patterns(patterns);

    // output pattern set
    if (savePatterns)
    {
        ofstream patOut("patterns.txt");
        for (auto &p : patterns)
        {
            for (auto &c : p)
            {
                patOut << c;
            }
            if (p != patterns.back())
            {
                patOut << endl;
            }
        }
    }

    /* -------------------- */
    /* parsing and calculating spacedWords */
    vector<Species> species;
    if (inFiles.empty() && !input_filename.empty())
    {
        cout << " --------------\nDetected Multifasta!\n";
        parser(input_filename, species);
    }
    else if (!inFiles.empty())
    {
        cout << " --------------\nDetected input folder!\n";
        sw_parser(inFiles, species, patterns);
    }

    /* -------------------- */
    /* calculating spaced-words */
    cout << "Calculating spaced-words...\n";
    double start_sw = wall_time();
    if (inFiles.empty() && !input_filename.empty())
    {
        for (unsigned int i = 0; i < species.size(); ++i)
        {
            for (const auto &pattern : patterns)
            {
                vector<Word> sw;
                spacedWords(species[i], pattern, sw);
            }
        }
    }
    else if (!inFiles.empty())
    {
        for (unsigned int i = 0; i < species.size(); ++i)
        {
            for (const auto &pattern : patterns)
            {
                spacedWords(species[i], pattern);
            }
        }
    }
    double end_sw = wall_time();
    cout << "Tiempo de spaced-words: " << (end_sw - start_sw) << " s\n";

    /* -------------------- */
    /* calculating matches */
    cout << " --------------\nCalculating matches...\n";
    double start_matches = wall_time();
    vector<vector<double>> distance(species.size(), vector<double>(species.size()));

    for (unsigned int i = 0; i < species.size(); ++i)
    {
        distance[i][i] = 0;

        for (auto j = species.size() - 1; j > i; --j)
        {
            double mismatch_rate = calc_matches(species[i], species[j], weight, dc, threshold, patterns, outputScores);
            distance[i][j] = calc_distance(mismatch_rate);
            if (mismatch_rate > 0.8541)
                tooDistant = true;
            distance[j][i] = distance[i][j];
        }
    }

    double end_matches = wall_time();
    cout << "Tiempo de matches: " << (end_matches - start_matches) << " s\n";

    /* -------------------- */
    /* output distance matrix */
    outputDistanceMatrix(species, output_filename, distance, tooDistant);

    double end_total = wall_time();
    cout << " --------------\nTotal run-time:\n";
    cout << (end_total - start) << " s\n";

    return 0;
}
