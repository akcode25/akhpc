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
*/
