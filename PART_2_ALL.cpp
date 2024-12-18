/* 10. Multiply two square matrices (1000,2000 or 3000 dimensions). Compare 
the performance of a sequential and parallel algorithm using open MP.*/

#include <iostream>
#include <vector>
#include <omp.h>
#include <ctime>
using namespace std;

void multiply(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int N, bool parallel) {
    #pragma omp parallel for collapse(2) if(parallel)  // Parallelize if "parallel" is true
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < N; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
int main() {
    int N;
    cout << "Enter matrix size (e.g., 1000, 2000, 3000): "; cin >> N;
    vector<vector<int>> A(N, vector<int>(N, 1)), B(N, vector<int>(N, 1)), C(N, vector<int>(N, 0));

    clock_t start = clock();
    multiply(A, B, C, N, false);  // Sequential multiplication
    clock_t end = clock();
    cout << "Time taken for sequential multiplication: " << double(end - start) / CLOCKS_PER_SEC << " seconds" << endl;

    start = clock();
    multiply(A, B, C, N, true);   // Parallel multiplication
    end = clock();
    cout << "Time taken for parallel multiplication (OpenMP): " << double(end - start) / CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}

/*
$ gedit ak10.cpp
$ g++ -fopenmp ak10.cpp
$ ./a.out

Enter matrix size (e.g., 1000, 2000, 3000): 1000
Time taken for sequential multiplication: 6.82819 seconds
Time taken for parallel multiplication (OpenMP): 0.62969 seconds
*/

_________________________________________________________________________________________

/*6. Assume you have n robots which pick mangoes in a farm. 
WAPT calculate the total number of mangoes picked by n robots parallely using MPI.*/

#include <iostream>
#include <mpi.h>
#include <ctime>
#include <cstdlib>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);			// Initialize MPI environment
    int world_size, world_rank;     // Get number of processes and rank of current process
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Initialize random seed and generate random number of mangoes
    srand(static_cast<unsigned>(time(0)) + world_rank); //ensures each process generates a different random no.
    int mangoes_picked = rand() % 101;  //[0, 100] rand no. => no. of mangoes picked by the robot
    cout << "Robot " << world_rank << " picked " << mangoes_picked << " mangoes." << endl;

    int total_mangoes = 0;
    MPI_Reduce(&mangoes_picked, &total_mangoes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) 
        cout << "Total mangoes picked by all robots: " << total_mangoes << endl;

    MPI_Finalize();
    return 0;
}

/*
$ gedit ak6.cpp
$ mpic++ ak6.cpp
$ mpirun -np 5 ./a.out

Robot 2 picked 25 mangoes.
Robot 4 picked 6 mangoes.
Robot 0 picked 15 mangoes.
Robot 1 picked 29 mangoes.
Robot 3 picked 49 mangoes.
Total mangoes picked by all robots: 124
*/

_________________________________________________________________________________________


/*7. Design a program that implements application of MPI Collective Communications.*/

#include <iostream>
#include <mpi.h>
#include <vector>
#include <cstdlib>  //rand() and srand()

using namespace std;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);  // Determines total no. of processes in communicator MPI_COMM_WORLD, Stores this value in world_size.
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  // Determines the rank (ID) of current process in the communicator, Stores this value in world_rank.

	int n = 10;
	vector<int> local_array(n);
	int local_sum = 0, total_sum = 0;

	srand(static_cast<unsigned>(time(0)) + world_rank);
	for (int i = 0; i < n; ++i) {
    	local_array[i] = rand() % 100;   // [0, 99]
    	local_sum += local_array[i];
	}
	cout << "Process " << world_rank << " local array: ";
	for (int i : local_array) cout << i << " ";
	cout << "\nProcess " << world_rank << " local sum: " << local_sum << endl;

	MPI_Reduce(&local_sum, &total_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&total_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//Broadcasts a variable (total_sum) from the root process to all other processes.
	//After this, all processes have the same value of total_sum

	if (world_rank == 0) {
    	double average = static_cast<double>(total_sum) / (n * world_size);
    	cout << "Total sum: " << total_sum << endl;
    	cout << "Average: " << average << endl;
	}
	MPI_Finalize();
	return 0;
}

/*
$ gedit akhpc7.cpp
$ mpic++ akhpc7.cpp
$ mpirun -np 5 ./a.out

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

_________________________________________________________________________________________


// 8. Implement Cartesian Virtual Topology in MPI.

#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int dims[2], coords[2], north, south, east, west, value;
    int periods[2] = {0, 0}; // No wrap-around in both dimensions
    MPI_Comm cart_comm; //communicator

    MPI_Dims_create(world_size, 2, dims);  //2D grid
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm); //Creates a new communicator with Cartesian topology
                 //Original communicator         //reorder

    MPI_Cart_coords(cart_comm, world_rank, 2, coords);  //Get process coordinates in the Cartesian grid
                                            //Array-store coordinates.
    MPI_Cart_shift(cart_comm, 0, 1, &north, &south); // Get neighbors in X direction of cartesian topology
    MPI_Cart_shift(cart_comm, 1, 1, &west, &east); // Get neighbors in Y direction of cartesian topology
    //(0-vertical,1-horizontal) //Dist of the shift (find immediate neighbors)   
    value = world_rank; // Value to send/receive
    cout << "Process " << world_rank << " at (" << coords[0] << ", " << coords[1] << ") has value: " << value << endl;

    // Send and receive data with neighbors (if exists)
    if (north != MPI_PROC_NULL) MPI_Send(&value, 1, MPI_INT, north, 0, cart_comm);
    if (south != MPI_PROC_NULL) MPI_Recv(&value, 1, MPI_INT, south, 0, cart_comm, MPI_STATUS_IGNORE); //no need details abt the rcv opern.
    if (west != MPI_PROC_NULL) MPI_Send(&value, 1, MPI_INT, west, 0, cart_comm);
    if (east != MPI_PROC_NULL) MPI_Recv(&value, 1, MPI_INT, east, 0, cart_comm, MPI_STATUS_IGNORE);

    MPI_Finalize();
    return 0;
}

/*
$ gedit ak8.cpp
$ mpic++ ak8.cpp
$ mpirun -np 4 ./a.out

Process 0 at (0, 0) has value: 0
Process 1 at (0, 1) has value: 1
Process 2 at (1, 0) has value: 2
Process 3 at (1, 1) has value: 3
*/

_________________________________________________________________________________________


/*9. Design a MPI program that uses blocking send/receive 
routines and non blocking send/receive routines.*/

#include <iostream>
#include <mpi.h>
#include <vector>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int TAG = 0, DATA_SIZE = 10; //message tag to identify communication

    vector<int> send_data(DATA_SIZE, world_rank);   // Initialize send_data with world_rank values
    vector<int> recv_data(DATA_SIZE);               // Initialize empty recv_data

    MPI_Request send_request, recv_request;         // Request variables for non-blocking communication

    if (world_rank == 0) {

        // Blocking send => The function does not return until the send operation is completed.
        cout << "Process 0 sending data (blocking): ";
        for (int i : send_data) cout << i << " "; //prints the values (all initialized to 0).
        cout << endl;
        MPI_Send(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD);
        cout << "Process 0: Blocking send completed." << endl;

        // Non-blocking send => send is initiated, prog can continue executing other code until MPI_Wait
        MPI_Isend(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD, &send_request);
        cout << "Process 0 non-blocking send initiated." << endl;
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);  // Wait for the non-blocking send to complete
        cout << "Process 0 non-blocking send completed." << endl;
    } 
    else if (world_rank == 1) {
	
        // Blocking receive => Waits until the data is fully received.
        MPI_Recv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "\nProcess 1 received data (blocking): ";
        for (int i : recv_data) cout << i << " ";
        cout << endl;
        cout << "Process 1: Blocking receive completed." << endl;

        // Non-blocking receive
        MPI_Irecv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, &recv_request);
        cout << "Process 1 non-blocking receive initiated." << endl;
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);  // Wait for the non-blocking receive to complete
        cout << "Process 1 non-blocking receive completed: ";
        for (int i : recv_data) cout << i << " ";
        cout << endl;
    }

    MPI_Finalize();
    return 0;
}


/*
$ gedit ak9.cpp
$ mpic++ ak9.cpp
$ mpirun -np 9 ./a.out

Process 0 sending data (blocking): 0 0 0 0 0 0 0 0 0 0 
Process 0: Blocking send completed.
Process 0 non-blocking send initiated.
Process 0 non-blocking send completed.
Process 1 received data (blocking): 0 0 0 0 0 0 0 0 0 0 
Process 1: Blocking receive completed.
Process 1 non-blocking receive initiated.
Process 1 non-blocking receive completed: 0 0 0 0 0 0 0 0 0 0 
*/


_________________________________________________________________________________________

// 8th compact

#include<iostream>
#include<mpi.h>
using namespace std;
int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int dims[2], coords[2], north, south, east, west, value;
    int periods[2] = {0,0};
    MPI_Comm cart_comm;

    MPI_Dims_create(world_size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
    MPI_Cart_coords(cart_comm, world_rank, 2, coords);
    MPI_Cart_shift(cart_comm, 0, 1, &north, &south);
    MPI_Cart_shift(cart_comm, 1, 1, &west, &east);

    value = world_rank;
    cout<<"process "<<world_rank<<" at ("<<coords[0]<<", "<<coords[1]<<" ) has value: "<<value;

    if (north != MPI_PROC_NULL) MPI_Send(&value, 1, MPI_INT, north, 0, cart_comm);
    if (south != MPI_PROC_NULL) MPI_Recv(&value, 1, MPI_INT, south, 0, cart_comm, MPI_STATUS_IGNORE);
    if (west != MPI_PROC_NULL) MPI_Send(&value, 1, MPI_INT, west, 0, cart_comm);
    if (east != MPI_PROC_NULL) MPI_Recv(&value, 1, MPI_INT, east, 0, cart_comm, MPI_STATUS_IGNORE);

    MPI_Finalize();
    return 0;
}

// mpic++ ak8.cpp
// mpirun -np 4 ./a.out

_________________________________________________________________________________________

//9th compact

#include<iostream>
#include<mpi.h>
#include<vector>
using namespace std;
int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    const int TAG = 0, DATA_SIZE = 10;
    vector<int> send_data(DATA_SIZE, world_rank);
    vector<int> recv_data(DATA_SIZE);
    MPI_Request send_request, recv_request;

    if (world_rank==0){
        cout<<"Process 0 send initiated (blocking): ";
        for (int i:send_data) cout<<i<<" ";
        MPI_Send(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD);
        cout<<"Process 0 send completed (blocking).";
        
        MPI_Isend(send_data.data(), DATA_SIZE, MPI_INT, 1, TAG, MPI_COMM_WORLD, &send_request);
        cout<<"Process 0 send initiated (non-blocking).";
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
        cout<<"Process 0 send completed (non-blocking).";
    }
    else if (world_rank==1){
        MPI_Recv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout<<"Process 1 receive initiated (blocking): ";
        for (int i:recv_data) cout<<i<<" ";
        cout<<"Process 1 recieve completed(blocking).";
        
        MPI_Irecv(recv_data.data(), DATA_SIZE, MPI_INT, 0, TAG, MPI_COMM_WORLD, &recv_request);
        cout<<"Process 1 receive initiated (non-blocking).";
        MPI_Wait(&recv_request, MPI_STSTUS_IGNORE);
        cout<<"Process 1 recieve completed (non-blocking): ";
        for (int i:recv_data) cout<<i<<" ";
    }
    MPI_Finalize();
    return 0;
}

// mpi++ ak9.cpp
// mpirun -np 9 ./a.out


