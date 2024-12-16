#include <iostream>
#include <vector>
#include <omp.h>  // OpenMP
#include <ctime> 

void multiply(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int N, bool parallel) {
    #pragma omp parallel for collapse(2) if(parallel)  // Parallelize if "parallel" is true
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < N; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main() {
    int N;
    cout << "Enter matrix size (e.g., 1000, 2000, 3000): "; cin >> N;

    vector<vector<int>> A(N, vector<int>(N, 1)), B(N, vector<int>(N, 1)), C(N, vector<int>(N, 0));

    clock_t start = clock();
    multiply(A, B, C, N, false);  // Sequential multiplication
    clock_t end = clock();
    cout << "Time taken for sequential multiplication: " << double(end - start) / CLOCKS_PER_SEC << " seconds" << endl;

    start = clock();
    multiply(A, B, C, N, true);   // Parallel multiplication
    end = clock();
    cout << "Time taken for parallel multiplication (OpenMP): " << double(end - start) / CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}

/*
Enter matrix size (e.g., 1000, 2000, 3000): 1000
Time taken for sequential multiplication: 6.82819 seconds
Time taken for parallel multiplication (OpenMP): 0.62969 seconds
*/
