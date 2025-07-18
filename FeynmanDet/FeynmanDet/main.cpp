//
//  main.cpp
//  Path_simulator
//  build e depois: main.exe
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
#include <iostream>
#include <stdlib.h>
#include "BaseFunctionalty/circuit.h"
#include "Simulators/simulator_RG.hpp"
#include "Simulators/simulator_PB.hpp"
#include "Simulators/simulator_all.hpp"
#include <time.h>
#include "States.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

int n_threads=-1;

int main(int argc, char *argv[]) {

    char fileName[256];
    int algorithm = 1, fs_int=0;

    if (argc>=2) { // get circuit name
        snprintf (fileName, 255, "circuits/circuit_%s.data", argv[1]);
    }
    else {
        snprintf (fileName, 255, "circuits/circuit_19.data");
    }
    // get final state (integer)
    if (argc>=3) {
        fs_int = atoi(argv[2]);
    }
    // get algorithm
    // 0 - all paths
    // 1 - RG paths
    // 2 - PB paths
    if (argc>=4) {
        algorithm = atoi(argv[3]);
    }
    // get num threads (OpenMP)
    if (argc>=5) {
        n_threads = atoi(argv[4]);
    }

#if defined(_OPENMP)
    fprintf(stderr, "OpenMP enabled: %d cores.\n", omp_get_num_procs());
#endif
    
    fprintf (stdout, "Circuit: %s\n", fileName);
    switch (algorithm) {
        case 0:
            fprintf (stdout, "ALL_PATHS\n");
            break;
        case 1:
            fprintf (stdout, "RG_PATHS\n");
            break;
        case 2:
            fprintf (stdout, "PB_PATHS\n");
            break;

        default:
            fprintf (stdout, "UNKNOWN ALGORITHM. Default = RG_PATHS\n");
            algorithm = 1;
            break;
    }
    // Call the read_circuit function
    TCircuit * circuit = read_circuit(fileName);
    // Check if the function succeeded
    if (circuit == NULL) {
        fprintf(stderr, "ERROR: unable to load circuit from file %s!\n", fileName);
        return 0;
    }
 
    
    const int L = circuit->size->num_layers;
    printf("L= %d \n", L);
    const int NQ = circuit->size->num_qubits;
    if (NQ >64) {
        fprintf (stderr, "ERROR: maximum 64 qubits allowed\n");
        return 0;
    }
    
    const StateT N = 1 << NQ;
    printf("NQ= %d \n", NQ);
    printf("N= %llu \n", N);
    
    //print_circuit(circuit);
    //for (fs=0 ; fs < N ; fs++) {
        clock_t start, end;
        start=clock();


        StateT init_state = 0;
    StateT final_state= (StateT) fs_int;
        float aR=0;
        float aI=0;
        
        switch (algorithm) {
            case 0:
                simulate_all_paths(circuit, init_state, final_state, aR, aI);
                break;
            case 1:
                simulate_RG_paths(circuit, init_state, final_state, aR, aI);
                break;
            case 2:
                simulate_PB_paths(circuit, init_state, final_state, aR, aI);
                break;

            default:
                break;
        }
        
        end=clock();
        double time_taken=double(end - start)/double(CLOCKS_PER_SEC);
        printf("Time taken is: %lf secs\n", time_taken);
    //}
    // Free the allocated memory
    free(circuit->size);
    free(circuit->layers);
    free(circuit);

    printf("\n");
    return 0;
}
