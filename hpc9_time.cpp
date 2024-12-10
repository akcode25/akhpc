// no need
// not tested

#include <iostream>
#include <mpi.h>
#include <vector>
#include <ctime>

using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int TAG = 0;
    const int DATA_SIZE = 10;
    vector<int> send_data(DATA_SIZE, world_rank);
    vector<int> recv_data(DATA_SIZE);

    clock_t start_time, end_time; 
    double blocking_time = 0.0; 
    double non_blocking_time = 0.0;  

    if (world_rank == 0) {
        cout << "Process 0 sending data: ";
        for (int i : send_data) cout << i << " ";
        cout << endl;
        
        start_time = clock();  
        MPI_Send(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD);
        end_time = clock(); 
        blocking_time += double(end_time - start_time) / CLOCKS_PER_SEC;
    } 
    else if (world_rank == 1) {
        start_time = clock(); 
        MPI_Recv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        end_time = clock(); 
        blocking_time += double(end_time - start_time) / CLOCKS_PER_SEC; 
        cout << "Process 1 received data: ";
        for (int i : recv_data) cout << i << " ";
        cout << endl;
    }

    MPI_Request send_request, recv_request;
    
    if (world_rank == 0) {
        start_time = clock();  
        MPI_Isend(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD, &send_request);
        cout << "Process 0 non-blocking send initiated." << endl;
        MPI_Wait(&send_request, MPI_STATUS_IGNORE); 
        end_time = clock(); 
        non_blocking_time += double(end_time - start_time) / CLOCKS_PER_SEC; 
    } 
    else if (world_rank == 1) {
        start_time = clock();
        MPI_Irecv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, &recv_request);
        cout << "Process 1 non-blocking receive initiated." << endl;
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE); 
        end_time = clock();  
        non_blocking_time += double(end_time - start_time) / CLOCKS_PER_SEC; 
        cout << "Process 1 non-blocking receive completed: ";
        for (int i : recv_data) cout << i << " ";
        cout << endl;
    }

    if (world_rank == 0) {
        cout << "Total time taken for blocking communication: " << blocking_time << " seconds" << endl;
        cout << "Total time taken for non-blocking communication: " << non_blocking_time << " seconds" << endl;
    }

    MPI_Finalize();
    return 0;
}

/*
//not tested:
Process 0 sending data: 0 0 0 0 0 0 0 0 0 0 
Process 1 received data: 0 0 0 0 0 0 0 0 0 0
Process 0 non-blocking send initiated.
Process 0 non-blocking send completed.
Process 1 non-blocking receive initiated.
Process 1 non-blocking receive completed: 0 0 0 0 0 0 0 0 0 0
Total time taken for blocking communication: 0.000234 seconds
Total time taken for non-blocking communication: 0.000189 seconds
*/
