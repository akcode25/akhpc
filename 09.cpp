#include <iostream>
#include <mpi.h>
#include <vector>

using namespace std;
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int TAG = 0, DATA_SIZE = 10;
    vector<int> send_data(DATA_SIZE, world_rank); 
    vector<int> recv_data(DATA_SIZE);              

    MPI_Request send_request, recv_request;

    if (world_rank == 0) {
        cout << "Process 0 sending data (blocking): ";
        for (int i : send_data) cout << i << " ";
        cout << endl;

        // Blocking send
        MPI_Send(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD);
        cout << "Process 0: Blocking send completed." << endl;

        // Non-blocking send
        MPI_Isend(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD, &send_request);
        cout << "Process 0 non-blocking send initiated." << endl;
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);  
        cout << "Process 0 non-blocking send completed." << endl;
    } 
    else if (world_rank == 1) {
			
        // Blocking receive
        MPI_Recv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "Process 1 received data (blocking): ";
        for (int i : recv_data) cout << i << " ";
        cout << endl;
        cout << "Process 1: Blocking receive completed." << endl;

        // Non-blocking receive
        MPI_Irecv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, &recv_request);
        cout << "Process 1 non-blocking receive initiated." << endl;
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE); 
        cout << "Process 1 non-blocking receive completed: ";
        for (int i : recv_data) cout << i << " ";
        cout << endl;
    }

    MPI_Finalize();
    return 0;
}

/*
mpic++ akhpc9.cpp
mpirun -np 9 ./a.out

Process 0 sending data (blocking): 0 0 0 0 0 0 0 0 0 0 
Process 0: Blocking send completed.
Process 0 non-blocking send initiated.
Process 0 non-blocking send completed.
Process 1 received data (blocking): 0 0 0 0 0 0 0 0 0 0 
Process 1: Blocking receive completed.
Process 1 non-blocking receive initiated.
Process 1 non-blocking receive completed: 0 0 0 0 0 0 0 0 0 0 
*/
