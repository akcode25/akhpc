#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

vector<int> fib;
bool ready = false;

void generate_fib(int limit) {
	fib.push_back(0), fib.push_back(1);
	for (int i = 2; i < limit; i++)
    	fib.push_back(fib[i - 1] + fib[i - 2]);
	ready = true;
}

void print_fib() {
	while (!ready);
	cout << "Fibonacci: ";
	for (int i : fib)
    	cout << i << " ";
	cout << endl;
}

int main() {
	int limit;
	cout << "Limit: ";
	cin >> limit;

	double start_serial = omp_get_wtime();
	fib.clear();
	generate_fib(limit);
	print_fib();
	double end_serial = omp_get_wtime();
	cout << "Serial Time: " << end_serial - start_serial << " seconds\n";

	double start_parallel = omp_get_wtime();
	#pragma omp parallel
	{
    	#pragma omp single
    	{
        	fib.clear();
        	generate_fib(limit);
        	print_fib();
    	}
	}
	double end_parallel = omp_get_wtime();
	cout << "Parallel Time: " << end_parallel - start_parallel << " seconds\n";
    return 0;
}

/*
Limit: 15
Fibonacci: 0 1 1 2 3 5 8 13 21 34 55 89 144 233 377
Serial Time: 4.7139e-05 seconds
Fibonacci: 0 1 1 2 3 5 8 13 21 34 55 89 144 233 377
Parallel Time: 0.00110327 seconds
*/