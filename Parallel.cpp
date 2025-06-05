// Rehan Tariq
// 22i-0965
// CS-6A


#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

#define GAP_PENALTY -2
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1

using namespace std;

struct Sequence 
{
    string id;
    string seq;
};

// Reads sequences from a FASTA file.
vector<Sequence> loadSequences(const string &filename) 
{
    vector<Sequence> sequences;
    ifstream file(filename);
    if (!file) 
    {
        cerr << "Error: Unable to open file " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    string line;
    string currentID, currentSeq;
    
    while (getline(file, line)) 
    {
        if (line.empty()) continue;
        if (line[0] == '>') 
        {
            if (!currentID.empty()) 
            {
                sequences.push_back({currentID, currentSeq});
                currentSeq.clear();
            }
            currentID = line.substr(1); // Remove '>' from header
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
    
    file.close();
    return sequences;
}

// Needleman-Wunsch alignment to calculate alignment score between two sequences.
int computeAlignmentScore(const string &s1, const string &s2) 
{
    int n = s1.size(), m = s2.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    
    for (int i = 0; i <= n; i++)
        dp[i][0] = i * GAP_PENALTY;
        
    for (int j = 0; j <= m; j++)
        dp[0][j] = j * GAP_PENALTY;
        
    for (int i = 1; i <= n; i++) 
    {
        for (int j = 1; j <= m; j++) 
        {
            int diag = dp[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_SCORE);
            int up = dp[i - 1][j] + GAP_PENALTY;
            int left = dp[i][j - 1] + GAP_PENALTY;
            dp[i][j] = max(diag, max(up, left));
        }
    }
    return dp[n][m];
}

// Distributes pairwise computations among MPI processes and gathers the results.
void computeDistributedDistanceMatrix(const vector<Sequence> &seqs, int rank, int size) 
{
    int totalSeqs = seqs.size();
    int totalPairs = (totalSeqs * (totalSeqs - 1)) / 2;
    int pairsPerProc = totalPairs / size;
    int remainder = totalPairs % size;

    vector<vector<double>> localMatrix(totalSeqs, vector<double>(totalSeqs, 0.0));
    double tStart = MPI_Wtime();

    int startPair = rank * pairsPerProc + min(rank, remainder);
    int endPair = startPair + pairsPerProc + (rank < remainder ? 1 : 0);

    int pairCount = 0;
    for (int i = 0; i < totalSeqs; i++) 
    {
        for (int j = i + 1; j < totalSeqs; j++) 
        {
            if (pairCount >= startPair && pairCount < endPair) 
            {
                int score = computeAlignmentScore(seqs[i].seq, seqs[j].seq);
                int maxLen = max(seqs[i].seq.size(), seqs[j].seq.size());
                double dist = 1.0 - (static_cast<double>(score) / maxLen);
                localMatrix[i][j] = dist;
                localMatrix[j][i] = dist; // Ensure symmetry
            }
            pairCount++;
        }
    }

    // Gather all local matrices into one global matrix.
    vector<double> globalMatrix(totalSeqs * totalSeqs, 0.0);
    for (int i = 0; i < totalSeqs; i++) 
    {
        MPI_Allreduce(localMatrix[i].data(), &globalMatrix[i * totalSeqs], totalSeqs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    double tEnd = MPI_Wtime();
    double elapsed = tEnd - tStart;

    if (rank == 0)
    {
        cout << "Distance Matrix:\n\t";
        for (const auto &s : seqs) cout << s.id << "\t";
        cout << "\n";
        
        for (int i = 0; i < totalSeqs; i++) 
        {
            cout << seqs[i].id << "\t";
            for (int j = 0; j < totalSeqs; j++) 
            {
                cout << fixed << setprecision(6) << globalMatrix[i * totalSeqs + j] << "\t";
            }
            cout << "\n";
        }
        
        cout << "\nExecution Time: " << fixed << setprecision(6) << elapsed << " seconds\n";
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    string fastaFile = "data2.txt";
    vector<Sequence> sequences = loadSequences(fastaFile);

    computeDistributedDistanceMatrix(sequences, rank, size);

    MPI_Finalize();
    return 0;
}

