#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

int main() {
	int num_students;
	cout << "Enter number of students: ";
	cin >> num_students;

	vector<double> CGPA(num_students);
	cout << "Enter the CGPAs of the students:\n";
	for (double &cgpa : CGPA) cin >> cgpa;

	double max_cgpa_serial = CGPA[0], max_cgpa_parallel = CGPA[0];

	// Serial execution
	double start_serial = omp_get_wtime();
	for (int i = 0; i < num_students; i++)
    	if (CGPA[i] > max_cgpa_serial) max_cgpa_serial = CGPA[i];
	double end_serial = omp_get_wtime();

	// Parallel execution
	double start_parallel = omp_get_wtime();
	#pragma omp parallel for shared(CGPA, max_cgpa_parallel)
	for (int i = 0; i < num_students; i++) {
    	#pragma omp critical
    	{
        	if (CGPA[i] > max_cgpa_parallel)
            	max_cgpa_parallel = CGPA[i];
    	}
	}
	double end_parallel = omp_get_wtime();

	cout << "Highest CGPA (Serial): " << max_cgpa_serial << "\n";
	cout << "Serial Time: " << end_serial - start_serial << " seconds\n";
	cout << "Highest CGPA (Parallel): " << max_cgpa_parallel << "\n";
	cout << "Parallel Time: " << end_parallel - start_parallel << " seconds\n";
	return 0;
}

/*
Enter number of students: 5
Enter the CGPAs of the students:
9.5
8.6
9.7
9.3
8.8
*/
