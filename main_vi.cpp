//
//  main_vi.cpp
//  Path_simulator
//build e depois: main_vi.exe
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
#include <iostream>
#include <stdlib.h>
#include "circuit.h"
#include <time.h>

#include "simulador.cpp"


   int main() {
    const char * fileName = "teste3.data"; 

    // Call the read_circuit function
    TCircuit * circuit = read_circuit(fileName);
    
    
    const int L = circuit->size->num_layers;
    printf("Lt= %d \n", L);
    const int NQ = circuit->size->num_qubits; 
    const long long N = (long long)(powf(2.f, (float)NQ));
    printf("NQ= %d \n", NQ);
    printf("N= %llu \n", N);
    
    clock_t start, end;
    start=clock();


    // Check if the function succeeded
    if (circuit != NULL) {
        //print_circuit(circuit);
        int* states = new int[(L-1)*NQ];
        int init_state;
        int final_state;
        for(int i=0;i<1; i++){
        	for(int f=0;f<N;f++){ //N
        		float aR=0;
				float aI=0;
				//simulate_all_paths(circuit, i, f, aR, aI);
                simulate_useful_paths(circuit, i, f, aR, aI);	
			}
		}

   		end=clock();
   		double time_taken=double(end - start)/double(CLOCKS_PER_SEC);
   		printf("Time taken is: %lf \n", time_taken);
        // Free the allocated memory
        free(circuit->size);
        free(circuit->layers);
        free(circuit);
    } else {
        fprintf(stderr, "Error reading circuit from file.\n");
    }

    printf("\n");
    system("pause");
    return 0;
}
