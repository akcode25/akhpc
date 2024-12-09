// 8. Implement Cartesian Virtual Topology in MPI.

#include <iostream>
#include <mpi.h>
#include <vector>
using namespace std;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int dims[2] = {0, 0};
	MPI_Dims_create(world_size, 2, dims);

	MPI_Comm cart_comm;
	int period[2] = {0, 0}; // Periodicity in both dimensions
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, true, &cart_comm);

	int coords[2];
	MPI_Comm_rank(cart_comm, &world_rank);
	MPI_Cart_coords(cart_comm, world_rank, 2, coords);

	int north, south, east, west;
	MPI_Cart_shift(cart_comm, 0, 1, &north, &south);
	MPI_Cart_shift(cart_comm, 1, 1, &west, &east);

	int value = world_rank;
	cout << "Process " << world_rank << " at (" << coords[0] << ", " << coords[1] << ") has value: " << value << endl;

	if (north != MPI_PROC_NULL) 
    	MPI_Send(&value, 1, MPI_INT, north, 0, cart_comm);
	if (south != MPI_PROC_NULL) 
    	MPI_Recv(&value, 1, MPI_INT, south, 0, cart_comm, MPI_STATUS_IGNORE);
	

	if (west != MPI_PROC_NULL) 
    	MPI_Send(&value, 1, MPI_INT, west, 0, cart_comm);
	if (east != MPI_PROC_NULL) 
    	MPI_Recv(&value, 1, MPI_INT, east, 0, cart_comm, MPI_STATUS_IGNORE);
	
	MPI_Finalize();
	return 0;
}

/*
//OUTPUT:

student@student-HP-Pro-Tower-280-G9-E-PCI-Desktop-PC:~$ gedit akhpc8.cpp
student@student-HP-Pro-Tower-280-G9-E-PCI-Desktop-PC:~$ mpic++ akhpc8.cpp
student@student-HP-Pro-Tower-280-G9-E-PCI-Desktop-PC:~$ mpirun -np 4 ./a.out
Invalid MIT-MAGIC-COOKIE-1 key
Process 0 at (0, 0) has value: 0
Process 1 at (0, 1) has value: 1
Process 2 at (1, 0) has value: 2
Process 3 at (1, 1) has value: 3
*/