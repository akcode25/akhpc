#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
	
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int dims[2], coords[2], north, south, east, west, value;
    int periods[2] = {0, 0}; 
    MPI_Comm cart_comm;

    MPI_Dims_create(world_size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
    
    MPI_Cart_coords(cart_comm, world_rank, 2, coords);

    MPI_Cart_shift(cart_comm, 0, 1, &north, &south);
    MPI_Cart_shift(cart_comm, 1, 1, &west, &east); 

    value = world_rank; 
    cout << "Process " << world_rank << " at (" << coords[0] << ", " << coords[1] << ") has value: " << value << endl;

    if (north != MPI_PROC_NULL) MPI_Send(&value, 1, MPI_INT, north, 0, cart_comm);
    if (south != MPI_PROC_NULL) MPI_Recv(&value, 1, MPI_INT, south, 0, cart_comm, MPI_STATUS_IGNORE);
    if (west != MPI_PROC_NULL) MPI_Send(&value, 1, MPI_INT, west, 0, cart_comm);
    if (east != MPI_PROC_NULL) MPI_Recv(&value, 1, MPI_INT, east, 0, cart_comm, MPI_STATUS_IGNORE);

    MPI_Finalize();
    return 0;
}

/*
gedit ak8.cpp
mpic++ ak8.cpp
mpirun -np 4 ./a.out

Process 0 at (0, 0) has value: 0
Process 1 at (0, 1) has value: 1
Process 2 at (1, 0) has value: 2
Process 3 at (1, 1) has value: 3
*/
