/*6. Assume you have n robots which pick mangoes in a farm. 
Wapt calculate the total number of mangoes picked by n robots parallely using MPI.*/

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);			// Initialize MPI environment

    // Get number of processes and rank of current process
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Initialize random seed and generate random number of mangoes
    srand(static_cast<unsigned>(time(0)) + world_rank);   //srand(world_rank);
    int mangoes_picked = rand() % 101;

    cout << "Robot " << world_rank << " picked " << mangoes_picked << " mangoes." << endl;

    // Use MPI_Reduce to calculate the total mangoes picked
    int total_mangoes = 0;
    MPI_Reduce(&mangoes_picked, &total_mangoes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Display the total mangoes picked by all robots (only on root process)
    if (world_rank == 0) {
        cout << "Total mangoes picked by all robots: " << total_mangoes << endl;
    }

    MPI_Finalize();
    return 0;
}




// #include <mpi.h>
// #include <iostream>
// #include <cstdlib>
// #include <ctime>

// using namespace std;

// int main(int argc, char** argv) {
// 	MPI_Init(&argc, &argv);        // Initialize MPI environment
// 	int world_size;
// 	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
// 	int world_rank;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
// 	srand(static_cast<unsigned>(time(0)) + world_rank);   //or, srand((unsigned int)time(0) + world_rank);
// 	int mangoes_picked = rand() % 101;
// 	cout << "Robot " << world_rank << " picked " << mangoes_picked << " mangoes." << endl;
// 	int total_mangoes;
// 	MPI_Reduce(&mangoes_picked, &total_mangoes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
// 	if (world_rank == 0) {
//     	cout << "Total mangoes picked by all robots: " << total_mangoes << endl;
// 	}
// 	MPI_Finalize();
// 	return 0;
// }

/*
//OUTPUT:
Robot 2 picked 25 mangoes.
Robot 4 picked 6 mangoes.
Robot 0 picked 15 mangoes.
Robot 1 picked 29 mangoes.
Robot 3 picked 49 mangoes.
Total mangoes picked by all robots: 124
*/