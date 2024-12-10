#include <iostream>
#include <mpi.h>
#include <vector>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int n = 10;
	vector<int> local_array(n);
	int local_sum = 0, total_sum = 0;

	srand(static_cast<unsigned>(time(0)) + world_rank);
	for (int i = 0; i < n; ++i) {
    	local_array[i] = rand() % 100;
    	local_sum += local_array[i];
	}
	cout << "Process " << world_rank << " local array: ";
	for (int i : local_array) cout << i << " ";
	cout << "\nProcess " << world_rank << " local sum: " << local_sum << endl;

	MPI_Reduce(&local_sum, &total_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&total_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (world_rank == 0) {
    	double average = static_cast<double>(total_sum) / (n * world_size);
    	cout << "Total sum: " << total_sum << endl;
    	cout << "Average: " << average << endl;
	}
	MPI_Finalize();
	return 0;
}

/*
//OUTPUT:

gedit akhpc7.cpp
mpic++ akhpc7.cpp
mpirun -np 5 ./a.out

Process 0 local array: 44 81 71 89 74 44 98 41 1 5
Process 0 local sum: 548
Process 3 local array: 62 42 67 68 15 31 40 14 96 2
Process 3 local sum: 437
Process 4 local array: 36 65 76 67 42 93 64 14 68 25
Process 4 local sum: 550
Process 2 local array: 36 8 70 74 74 29 58 11 69 79
Process 2 local sum: 508
Process 1 local array: 13 99 77 4 54 82 64 48 5 52
Process 1 local sum: 498
Total sum: 2541
Average: 50.82
*/