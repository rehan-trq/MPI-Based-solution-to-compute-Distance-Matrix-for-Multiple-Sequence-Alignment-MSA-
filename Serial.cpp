//Rehan Tariq
//22i-0965
//CS-6A 


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <chrono>

#define GAP_PENALTY -2
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1

using namespace std;
using namespace chrono;

struct Sequence 
{
    string id;
    string seq;
};

// Reads sequences from a FASTA file.
vector<Sequence> readFastaFile(const string &filename) 
{
    vector<Sequence> sequences;
    ifstream infile(filename);
    
    if (!infile) 
    {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }
    string line, currentID, currentSeq;
    
    while(getline(infile, line)) 
    {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (!currentID.empty()) 
            {
                sequences.push_back({currentID, currentSeq});
                currentSeq.clear();
            }
            currentID = line.substr(1);
        } 
        
        else 
        {
            currentSeq += line;
        }
    }
    
    if (!currentID.empty()) 
    {
        sequences.push_back({currentID, currentSeq});
    }
    
    infile.close();
    return sequences;
}

// Needleman-Wunsch algorithm to compute the alignment score between two sequences.
int needlemanWunsch(const string &s1, const string &s2) 
{
    int n = s1.size(), m = s2.size();
    vector<vector<int>> dp(n+1, vector<int>(m+1, 0));
    
    for (int i = 0; i <= n; i++)
        dp[i][0] = i * GAP_PENALTY;
        
    for (int j = 0; j <= m; j++)
        dp[0][j] = j * GAP_PENALTY;
        
    for (int i = 1; i <= n; i++) 
    {
        for (int j = 1; j <= m; j++) 
        {
            int diag = dp[i-1][j-1] + (s1[i-1] == s2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
            int up = dp[i-1][j] + GAP_PENALTY;
            int left = dp[i][j-1] + GAP_PENALTY;
            dp[i][j] = max(diag, max(up, left));
        }
    }
    return dp[n][m];
}

// Normalizes the raw score into a distance in the range 0.0 to 1.0.
double normalizeScore(int score, int minScore, int maxScore) 
{
    return (score - minScore) / static_cast<double>(maxScore - minScore);
}

// Computes and outputs the distance matrix.
void computeDistanceMatrix(const vector<Sequence> &sequences) 
{
    int numSeqs = sequences.size();
    vector<vector<int>> rawScores(numSeqs, vector<int>(numSeqs, 0));
    int minScore = numeric_limits<int>::max();
    int maxScore = numeric_limits<int>::min();

    auto tStart = high_resolution_clock::now();

    for (int i = 0; i < numSeqs; i++) 
    {
        for (int j = i + 1; j < numSeqs; j++) 
        {
            int score = needlemanWunsch(sequences[i].seq, sequences[j].seq);
            rawScores[i][j] = score;
            rawScores[j][i] = score;
            minScore = min(minScore, score);
            maxScore = max(maxScore, score);
        }
    }
    vector<vector<double>> distMatrix(numSeqs, vector<double>(numSeqs, 0.0));
    
    for (int i = 0; i < numSeqs; i++)
    {
        for (int j = 0; j < numSeqs; j++) 
        {
            if (i != j)
                distMatrix[i][j] = normalizeScore(rawScores[i][j], minScore, maxScore);
        }
    }

    auto tEnd = high_resolution_clock::now();
    double elapsed = duration<double>(tEnd - tStart).count();
    cout << "Distance Matrix:\n\t";
    
    for (const auto &s : sequences)
        cout << s.id << "\t";
    cout << "\n";

    for (int i = 0; i < numSeqs; i++) 
    {
        cout << sequences[i].id << "\t";
        for (int j = 0; j < numSeqs; j++) 
        {
            cout << fixed << setprecision(3) << distMatrix[i][j] << "\t";
        }
        cout << "\n";
    }
    cout << "\nExecution Time: " << fixed << setprecision(6) << elapsed << " seconds\n";
    ofstream outfile("distance_matrix_serial.txt");
    outfile << "Distance Matrix:\n\t";
    
    for (const auto &s : sequences)
        outfile << s.id << "\t";
        
    outfile << "\n";
    for (int i = 0; i < numSeqs; i++) 
    {
        outfile << sequences[i].id << "\t";
        for (int j = 0; j < numSeqs; j++) 
        {
            outfile << fixed << setprecision(3) << distMatrix[i][j] << "\t";
        }
        outfile << "\n";
    }
    
    outfile << "\nExecution Time: " << fixed << setprecision(6) << elapsed << " seconds\n";
    outfile.close();
}

int main(int argc, char **argv) 
{
    
    string filename = "data2.txt";
    vector<Sequence> sequences = readFastaFile(filename);
    
    if (sequences.empty()) 
    {
        cerr << "Error: No sequences found in " << filename << endl;
        return 1;
    }
    
    computeDistanceMatrix(sequences);
    return 0;
}


