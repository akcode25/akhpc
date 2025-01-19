g++ -fopenmp prog.cpp

/*1. */

#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

int main() {
	int n;
	cout << "Enter matrix size: ";
	cin >> n;
	vector<vector<int>> A(n, vector<int>(n));
	vector<int> x(n), y_serial(n, 0), y_parallel(n, 0);

	cout << "Enter matrix A :\n";
	for (auto &row : A) for (int &a : row) cin >> a;
	cout << "Enter vector x: \n";
	for (int &xi : x) cin >> xi;

	double start_serial = omp_get_wtime();
	for (int i = 0; i < n; i++)
    	    for (int j = 0; j < n; j++)
        	y_serial[i] += A[i][j] * x[j];
	double end_serial = omp_get_wtime();

	double start_parallel = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < n; i++)
    	    for (int j = 0; j < n; j++)
        	y_parallel[i] += A[i][j] * x[j];
	double end_parallel = omp_get_wtime();

	cout << "Serial Result y = ";
	for (int yi : y_serial) cout << yi << " ";
	cout << "\nSerial Time: " << end_serial - start_serial << " seconds\n";

	cout << "Parallel Result y = ";
	for (int yi : y_parallel) cout << yi << " ";
	cout << "\nParallel Time: " << end_parallel - start_parallel << " seconds\n";
    
	return 0;
}

/*
Enter matrix size: 3
Enter matrix A :
1 2 3
4 5 6
7 8 9
Enter vector x:
1 2 3

Serial Result y = 14 32 50
Serial Time: 9.91e-07 seconds
Parallel Result y = 14 32 50
Parallel Time: 0.00141599 seconds
*/

/*2. */

#include <iostream>
#include <vector>
#include <omp.h>
using namespace std;

int main() {
	vector<string> sections = {"Clothing", "Gaming", "Grocery", "Stationery"};
	vector<int> prices_parallel(sections.size(), 0), prices_serial(sections.size(), 0);

	// Serial execution
	cout<<"Serial execution: \n";
	double start_serial = omp_get_wtime();
	for (int i = 0; i < sections.size(); ++i) {
    	int num_items, total = 0;
    	cout << "Enter items & prices for " << sections[i] << " (Serial):\n";
    	cin >> num_items;
    	for (int j = 0; j < num_items; ++j) {
        	int price;
        	cin >> price;
        	total += price;
    	}
    	prices_serial[i] = total;
	}
	double end_serial = omp_get_wtime();

	// Parallel execution
	cout<<"\nParallel execution:\n";
	double start_parallel = omp_get_wtime();
	for (int i = 0; i < sections.size(); ++i) {
    	int num_items, total = 0;
    	cout << "Enter items & prices for " << sections[i] << " (Parallel):\n";
    	cin >> num_items;

    	#pragma omp parallel for reduction(+:total)
    	for (int j = 0; j < num_items; ++j) {
        	int price;
        	cin >> price;
        	total += price;
    	}
    	prices_parallel[i] = total;
	}
	double end_parallel = omp_get_wtime();

	// Final summary
	cout << "\nSerial Prices:\n";
	int overall_serial = 0;
	for (int i = 0; i < sections.size(); ++i) {
    	cout << sections[i] << ": " << prices_serial[i] << "\n";
    	overall_serial += prices_serial[i];
	}
	cout << "Overall Cost (Serial): " << overall_serial << "\n";
	cout << "Serial Time: " << end_serial - start_serial << " seconds\n";

	cout << "\nParallel Prices:\n";
	int overall_parallel = 0;
	for (int i = 0; i < sections.size(); ++i) {
    	cout << sections[i] << ": " << prices_parallel[i] << "\n";
    	overall_parallel += prices_parallel[i];
	}
	cout << "Overall Cost (Parallel): " << overall_parallel << "\n";
	cout << "Parallel Time: " << end_parallel - start_parallel << " seconds\n";

	return 0;
}


/*
Serial execution:
Enter items & prices for Clothing (Serial):
3
500 600 700
Enter items & prices for Gaming (Serial):
2
1500 2500
Enter items & prices for Grocery (Serial):
4
100 200 300 400
Enter items & prices for Stationery (Serial):
3
50 60 70

Parallel execution:
Enter items & prices for Clothing (Parallel):
3
500 600 700
Enter items & prices for Gaming (Parallel):
2
1500 2500
Enter items & prices for Grocery (Parallel):
4
100 200 300 400
Enter items & prices for Stationery (Parallel):
3
50 60 70

Serial Prices:
Clothing: 1800
Gaming: 4000
Grocery: 1000
Stationery: 180
Overall Cost (Serial): 6980
Serial Time: 19.3207 seconds

Parallel Prices:
Clothing: 1800
Gaming: 4000
Grocery: 820
Stationery: 180
Overall Cost (Parallel): 6800
Parallel Time: 19.8164 seconds
*/

/*3. */

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

/*4. */

#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

vector<int> fib;
bool ready = false;

void generate_fib(int limit) {
	fib.push_back(0), fib.push_back(1);
	for (int i = 2; i < limit; i++)
    	fib.push_back(fib[i - 1] + fib[i - 2]);
	ready = true;
}

void print_fib() {
	while (!ready);  //Busy-waits until ready becomes true, synchronization bw generation & printing processes.
	cout << "Fibonacci: ";
	for (int i : fib) cout << i << " "; cout << endl;
}

int main() {
	int limit;
	cout << "Limit: "; cin >> limit;

	double start_serial = omp_get_wtime();
	fib.clear();
	generate_fib(limit);
	print_fib();
	double end_serial = omp_get_wtime();
	cout << "Serial Time: " << end_serial - start_serial << " seconds\n";

	double start_parallel = omp_get_wtime();
	#pragma omp parallel
	{
    	#pragma omp single  //Ensures only a single thread executes this code block, while others wait.
    	{
        	fib.clear();
        	generate_fib(limit);
        	print_fib();
    	}
	}
	double end_parallel = omp_get_wtime();
	cout << "Parallel Time: " << end_parallel - start_parallel << " seconds\n";
    
	return 0;
}

/*
Limit: 15
Fibonacci: 0 1 1 2 3 5 8 13 21 34 55 89 144 233 377
Serial Time: 4.7139e-05 seconds
Fibonacci: 0 1 1 2 3 5 8 13 21 34 55 89 144 233 377
Parallel Time: 0.00110327 seconds
*/

/*5 */

/* 5. University awards gold medals to the student who has scored highest CGPA. WAPT TF the student with highest CGPA in a list of numbers using OpenMP.
*/

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

	double start_serial = omp_get_wtime();
	for (int i = 0; i < num_students; i++)
    	if (CGPA[i] > max_cgpa_serial) max_cgpa_serial = CGPA[i];
	double end_serial = omp_get_wtime();

	double start_parallel = omp_get_wtime();
	
	#pragma omp parallel for shared(CGPA, max_cgpa_parallel)
	for (int i = 0; i < num_students; i++) {
    	#pragma omp critical
    	{
        	if (CGPA[i] > max_cgpa_parallel) max_cgpa_parallel = CGPA[i];
    	}
	}

	// #pragma omp parallel for reduction(max: max_cgpa_p)
        // for (int i=0; i<n; i++)
        //     if(CGPA[i]>max_cgpa_p) max_cgpa_p=CGPA[i];
	
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
Highest CGPA (Serial): 9.7
Serial Time: 7.67e-07 seconds
Highest CGPA (Parallel): 9.7
Parallel Time: 0.000998799 seconds
*/

//_______________________________________________________________________________________________________

/*10. */

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
    multiply(A, B, C, N, false);  
    clock_t end = clock();
    cout << "Time taken for sequential multiplication: " << double(end - start) / CLOCKS_PER_SEC << " seconds" << endl;

    start = clock();
    multiply(A, B, C, N, true);  
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

/*6. */

#include <iostream>
#include <mpi.h>
#include <ctime>
#include <cstdlib>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);			
    int world_size, world_rank;    
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(static_cast<unsigned>(time(0)) + world_rank); 
    int mangoes_picked = rand() % 101;  //[0, 100] 
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

/*7. */

#include <iostream>
#include <mpi.h>
#include <vector>
#include <cstdlib>  //rand() and srand()

using namespace std;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);  // total no. of processes in communicator MPI_COMM_WORLD
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  // rank (ID) of current process in the communicator

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

/* 8. */

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

/*9. */

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


