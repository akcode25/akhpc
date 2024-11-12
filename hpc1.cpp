#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

int main() {
	int n;
	cout << "Enter matrix size: ";
	cin >> n;
	vector<vector<int>> A(n, vector<int>(n));
	vector<int> x(n), y_serial(n, 0), y_parallel(n, 0);

	cout << "Enter matrix A :\n";
	for (auto &row : A) for (int &a : row) cin >> a;
	cout << "Enter vector x: \n";
	for (int &xi : x) cin >> xi;

	// Serial execution
	double start_serial = omp_get_wtime();
	for (int i = 0; i < n; i++)
    	for (int j = 0; j < n; j++)
        	y_serial[i] += A[i][j] * x[j];
	double end_serial = omp_get_wtime();

	// Parallel execution
	double start_parallel = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < n; i++)
    	for (int j = 0; j < n; j++)
        	y_parallel[i] += A[i][j] * x[j];
	double end_parallel = omp_get_wtime();

	cout << "Serial Result y = ";
	for (int yi : y_serial) cout << yi << " ";
	cout << "\nSerial Time: " << end_serial - start_serial << " seconds\n";

	cout << "Parallel Result y = ";
	for (int yi : y_parallel) cout << yi << " ";
	cout << "\nParallel Time: " << end_parallel - start_parallel << " seconds\n";
    
	return 0;
}

/*
Enter matrix size: 3
Enter matrix A :
1 2 3
4 5 6
7 8 9
Enter vector x:
1 2 3
*/