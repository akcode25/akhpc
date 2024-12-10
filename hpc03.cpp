#include <iostream>
#include <omp.h>
using namespace std;

int main() {
	long long steps = 1000000000;  //10^9
	double step = 1.0 / steps, pi_serial = 0.0, pi_parallel = 0.0;

	// Serial execution
	double start_serial = omp_get_wtime();
	for (long long i = 0; i < steps; i++) {
    	double x = (i + 0.5) * step;
    	pi_serial += 4.0 / (1.0 + x * x);
	}
	pi_serial *= step;
	double end_serial = omp_get_wtime();

	// Parallel execution
	double start_parallel = omp_get_wtime();
	#pragma omp parallel
	{
    	double sum = 0.0;
    	#pragma omp for
    	for (long long i = 0; i < steps; i++) {
        	double x = (i + 0.5) * step;
        	sum += 4.0 / (1.0 + x * x);
    	}
    	#pragma omp critical
    	pi_parallel += sum * step;
	}
	double end_parallel = omp_get_wtime();

	cout << "Pi (Serial): " << pi_serial << "\n";
	cout << "Serial Time: " << end_serial - start_serial << " seconds\n";
	cout << "Pi (Parallel): " << pi_parallel << "\n";
	cout << "Parallel Time: " << end_parallel - start_parallel << " seconds\n";
	return 0;
}

/*
Pi (Serial): 3.14159
Serial Time: 5.66401 seconds
Pi (Parallel): 3.14159
Parallel Time: 0.325438 seconds
*/
