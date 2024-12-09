/* 10. Multiply two square matrices (1000,2000 or 3000 dimensions). Compare 
the performance of a sequential and parallel algorithm using open MP.*/

#include <iostream>
#include <vector>
#include <omp.h>  // OpenMP
#include <ctime> 

using namespace std;

// Function to multiply two matrices sequentially
void multiplySequential(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int N) {
	for (int i = 0; i < N; ++i) {
    	for (int j = 0; j < N; ++j) {
        	C[i][j] = 0;
        	for (int k = 0; k < N; ++k) {
            	C[i][j] += A[i][k] * B[k][j];
        	}
    	}
	}
}

// Function to multiply two matrices using OpenMP
void multiplyParallel(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int N) {
	#pragma omp parallel for collapse(2)
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
	int N; // Matrix size (1000, 2000, 3000)
	cout << "Enter matrix size (e.g., 1000, 2000, 3000): ";
	cin >> N;
    clock_t start, end;

	// Initialize matrices A, B, and C
	vector<vector<int>> A(N, vector<int>(N, 1)); // Matrix A with all elements = 1
	vector<vector<int>> B(N, vector<int>(N, 1)); // Matrix B with all elements = 1
	vector<vector<int>> C(N, vector<int>(N, 0)); // Matrix C to store result

	// Sequential multiplication
    start = clock();          // double start = omp_get_wtime();
	multiplySequential(A, B, C, N); 
    end = clock();            // double end = omp_get_wtime();
    double durationSeq = double(end - start) / CLOCKS_PER_SEC; 	// double durationSeq = end - start;
	cout << "Time taken for sequential multiplication: " << durationSeq << " seconds" << endl;

	// Parallel multiplication using OpenMP
	start = clock();                  //start = omp_get_wtime();
	multiplyParallel(A, B, C, N);
	end = clock();                    //end = omp_get_wtime();
	double durationPar = double(end - start) / CLOCKS_PER_SEC; //double durationPar = end - start;
	cout << "Time taken for parallel multiplication (OpenMP): " << durationPar << " seconds" << endl;

	return 0;
}

/*
//OUTPUT:
Enter matrix size (e.g., 1000, 2000, 3000): 1000
Time taken for sequential multiplication: 6.82819 seconds
Time taken for parallel multiplication (OpenMP): 0.62969 seconds
*/