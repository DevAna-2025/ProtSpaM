#include "calc_matches.h"
#include <cstdint>

using namespace std;

const int rows_col = 24;
int blosum62 [rows_col][rows_col] =
    {
        {4,  -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
        {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
        {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
        {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
        {0,  -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
        {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
        {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
        {0,  -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
        {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
        {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
        {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
        {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
        {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
        {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
        {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
        {1,  -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
        {0,  -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
        {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
        {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
        {0,  -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
        {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
        {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
        {0,  -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
        {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
    };

int multiMatch(const vector<Word> &sortedWords, int start) {
    int multiMatch_length = 1;
    while (sortedWords[start].key == sortedWords[start+1].key) {
        ++multiMatch_length;
        ++start;
    }
    return multiMatch_length;
}

// Posiciones don't-care del patron actual. Las calculo una vez por patron para
// no recorrer la mascara completa en cada comparacion de spaced-words.
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

// Calcula score y mismatches usando solo las posiciones don't-care ya filtradas.
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

// Precalculo de bloques con la misma key. Mantiene la misma idea de multiMatch,
// pero evita repetir esa busqueda dentro de los bucles principales.
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

void scoreOutput(map<int,int> &scores, string header1, string header2) {
    header1.erase(header1.find_last_not_of(" \n\r\t") + 1);
    header2.erase(header2.find_last_not_of(" \n\r\t") + 1);
    string fileName = "scores/" + header1 + "_" + header2 + ".scr";

    vector<pair<int, int>> scoreVector(scores.size());
    int i = 0;
    for (const auto &pair : scores) {
        scoreVector[i++] = pair;
    }
    sort(scoreVector.begin(), scoreVector.end());

    ofstream scoresOut;
    scoresOut.open(fileName);
    scoresOut << "score,freq\n";
    for (auto &s : scoreVector) {
        if (s.first > INT32_MIN) {
            scoresOut << s.first << "," << s.second << "\n";
        }
    }
}

double calc_mismatchRate (vector<pair<int,int> >& mismatchesDC) {
    int mismatches = 0;
    int dontCare = 0;
    for (auto pair : mismatchesDC) {
        mismatches += pair.first;
        dontCare += pair.second;
    }
    return (double) mismatches / dontCare;
}

double calc_matches (const Species& species1, const Species& species2, const int& weight, const int& dc,
                     const int& threshold, const vector<vector<char> >& patterns, const bool& outputScores) {
    unsigned int skip = 0;
    int total_mismatches = 0;
    int total_dc = 0;
    int score;
    int mismatches;
    map<int, int> scores;
    bool multi_done = false;
    const vector<char>& sequence1 = species1.seq;
    const vector<char>& sequence2 = species2.seq;
    vector<pair<int, int>> mismatchDontCare;

    for (unsigned int currentPattern = 0; currentPattern < patterns.size(); ++currentPattern) {
        const vector<Word>& spacedWords1 = species1.sorted_words[currentPattern];
        const vector<Word>& spacedWords2 = species2.sorted_words[currentPattern];
        const vector<char>& pattern = patterns[currentPattern];
        const vector<int> positions = dontcare_positions(pattern);
        const vector<int> run_lengths1 = key_run_lengths(spacedWords1);
        const vector<int> run_lengths2 = key_run_lengths(spacedWords2);
        for (unsigned int i = 0; i < spacedWords1.size(); ++i) {
            int bl1 = run_lengths1[i];
            if (bl1 > 1) {
                int best_score = threshold - 1;
                int best_mismatches = 0;
                unsigned int limit1 = i + bl1 - 1;
                for (; i <= limit1; ++i) {
                    bool singleMatch = true;
                    for (unsigned int j = skip; j < spacedWords2.size() && singleMatch; ++j) {
                        int bl2 = run_lengths2[j];
                        if (bl2 > 1) {
                            singleMatch = false;
                            unsigned int limit2 = j + bl2 - 1;
                            for (; j <= limit2; ++j) {
                                score = 0;
                                mismatches = 0;
                                if (spacedWords1[i].key > spacedWords2[j].key) {
                                    skip += bl2;
                                    continue;
                                }
                                if (spacedWords1[i].key < spacedWords2[j].key) {
                                    break;
                                }
                                if (spacedWords1[i].key == spacedWords2[j].key) {
                                    score_dontcare_positions(sequence1, sequence2,
                                                             spacedWords1[i].pos,
                                                             spacedWords2[j].pos,
                                                             positions, score, mismatches);
                                }
                                if (score >= threshold && score > best_score) {
                                    best_score = score;
                                    best_mismatches = mismatches;
                                }
                                if (i == limit1 && j == limit2) {
                                    skip += bl2;
                                    multi_done = true;
                                    break;
                                }
                            }
                        }
                        else {
                            if (spacedWords1[i].key > spacedWords2[j].key) {
                                skip = j + 1;
                                continue;
                            }
                            if (spacedWords1[i].key < spacedWords2[j].key) {
                                break;
                            }
                            if (spacedWords1[i].key == spacedWords2[j].key) {
                                score = 0;
                                mismatches = 0;
                                score_dontcare_positions(sequence1, sequence2,
                                                         spacedWords1[i].pos,
                                                         spacedWords2[j].pos,
                                                         positions, score, mismatches);
                                if (score >= threshold && score > best_score) {
                                    best_score = score;
                                    best_mismatches = mismatches;
                                }
                                if (i == limit1) {
                                    skip = j + 1;
                                    multi_done = true;
                                    break;
                                }
                                break;
                            }
                        }
                    }
                }
                if (outputScores) {
                    auto exists = scores.find(best_score);
                    if (exists != scores.end() ) {
                        scores[best_score] += 1;
                    }
                    else if (best_score){
                        scores.insert(pair<int, int>(best_score, 1));
                    }
                }
                if (best_score >= threshold) {
                    total_mismatches += best_mismatches;
                    total_dc += dc;
                }
                if (bl1 > 1 && multi_done) {
                    --i;
                    multi_done = false;
                }
            }
            else {
                bool go_on = true;
                for (unsigned int j = skip; j < spacedWords2.size() && go_on; ++j) {
                    int bl2 = run_lengths2[j];
                    if (bl2 > 1) {
                        go_on = false;
                        int best_score = threshold - 1;
                        int best_mismatches = 0;
                        unsigned int limit = j + bl2 - 1;
                        for (; j <= limit; ++j) {
                            score = 0;
                            mismatches = 0;
                            if (spacedWords1[i].key > spacedWords2[j].key) {
                                skip += bl2;
                                continue;
                            }
                            if (spacedWords1[i].key < spacedWords2[j].key) {
                                break;
                            }
                            if (spacedWords1[i].key == spacedWords2[j].key) {
                                score_dontcare_positions(sequence1, sequence2,
                                                         spacedWords1[i].pos,
                                                         spacedWords2[j].pos,
                                                         positions, score, mismatches);
                                if (score >= threshold && score > best_score) {
                                    best_score = score;
                                    best_mismatches = mismatches;
                                }
                                if (j == limit) {
                                    if (outputScores) {
                                        auto exists = scores.find(best_score);
                                        if (exists != scores.end() ) {
                                            scores[best_score] += 1;
                                        }
                                        else {
                                            scores.insert(pair<int, int>(best_score, 1));
                                        }
                                    }
                                    if (best_score > threshold) {
                                        total_mismatches += best_mismatches;
                                        total_dc += dc;
                                    }
                                    skip += bl2;
                                }
                            }
                        }
                    }
                    else {
                        if (spacedWords1[i].key > spacedWords2[j].key) {
                            skip = j + 1;
                            continue;
                        }
                        if (spacedWords1[i].key < spacedWords2[j].key) {
                            break;
                        }
                        if (spacedWords1[i].key == spacedWords2[j].key) {
                            skip = j + 1;
                            score = 0;
                            mismatches = 0;
                            score_dontcare_positions(sequence1, sequence2,
                                                     spacedWords1[i].pos,
                                                     spacedWords2[j].pos,
                                                     positions, score, mismatches);
                            if (outputScores) {
                                auto exists = scores.find(score);
                                if (exists != scores.end() ) {
                                    scores[score] += 1;
                                }
                                else {
                                    scores.insert(pair<int, int>(score, 1));
                                }
                            }
                            if (score >= threshold) {
                                total_mismatches += mismatches;
                                total_dc += dc;
                            }
                            break;
                        }
                    }
                }
            }
        }
        mismatchDontCare.emplace_back(total_mismatches, total_dc);
    }

    if (outputScores) {
        scoreOutput(scores, species1.header, species2.header);
    }

    return calc_mismatchRate(mismatchDontCare);
}
