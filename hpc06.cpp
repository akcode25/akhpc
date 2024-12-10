#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);	

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(static_cast<unsigned>(time(0)) + world_rank); //srand(world_rank); //or, srand((unsigned int)time(0) + world_rank);
    int mangoes_picked = rand() % 101;  
    cout << "Robot " << world_rank << " picked " << mangoes_picked << " mangoes." << endl;

    int total_mangoes = 0;
    MPI_Reduce(&mangoes_picked, &total_mangoes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0) cout << "Total mangoes picked by all robots: " << total_mangoes << endl;

    MPI_Finalize();
    return 0;
}

/*
Robot 2 picked 25 mangoes.
Robot 4 picked 6 mangoes.
Robot 0 picked 15 mangoes.
Robot 1 picked 29 mangoes.
Robot 3 picked 49 mangoes.
Total mangoes picked by all robots: 124
*/
