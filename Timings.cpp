#include <iostream>
#include <stdlib.h>

#include "circuit.h"
#include <math.h>
#include "complex.h"
#include "layer.hpp"
#include "layer.cpp"
#include "simulation_attempt.cpp"
#include "simulate_all.cpp"
#include "simulate_useful.cpp"
#include "circuit.cpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

int main() {
    const char* circuitNames[5] = {"circuit_1_size_3_depth_3.data", "circuit_2_size_4_depth_4.data", "circuit_3_size_3_depth_3.data","circuit_21_size_6_depth_5.data","circuit_211_size_7_depth_5.data"};
    const int numCircuits = 5;

    float durationValues[numCircuits];

    for (int circuitIndex = 0; circuitIndex < numCircuits; ++circuitIndex) {
        const char* fileName = circuitNames[circuitIndex];
		printf(fileName);
        // Call the read_circuit function
        TCircuit* circuit = read_circuit(fileName);
		
        const int L = circuit->size->num_layers;
        const int NQ = circuit->size->num_qubits;
        const int N = (int)(powf(2.f, (float)NQ));
        float avec[N];

        // Check if the function succeeded
        if (circuit != NULL) {
            // print_circuit(circuit);
			
            auto start_time = std::chrono::high_resolution_clock::now();
			
            int* states = new int[(L - 1) * NQ];
            int init_state;
            int final_state;
            for (int f = 0; f < N; f++) {
                float aR = 0;
                float aI = 0;
                avec[f] = simulate_useful_paths(circuit, 0, f, aR, aI);
            }
            
            std::ofstream outputFile("avec_output.txt");
            if (outputFile.is_open()) {
                for (int i = 0; i < N; i++) {
                    outputFile << avec[i] << "\n";
                }
                outputFile.close();
            } else {
                std::cerr << "Error opening output file.\n";
            }
			
            // Free the allocated memory
            free(circuit->size);
    		free(circuit->layers);
    		free(circuit);
			
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
            durationValues[circuitIndex] = duration;
            std::cout << durationValues[circuitIndex] << "\n";
        } else {
            fprintf(stderr, "Error reading circuit from file: %s\n", fileName);
            free(circuit->size);
    		free(circuit->layers);
    		free(circuit);
        }
    }
    

    // Output the duration values
    std::cout << "\nDuration values:\n";
    for (int i = 0; i < numCircuits; ++i) {
        std::cout << "Circuit " << i + 1 << ": " << durationValues[i] << " seconds\n";
    }
    
    std::ofstream outputFile("timing_output.txt");
            if (outputFile.is_open()) {
                for (int i = 0; i < numCircuits; i++) {
                    outputFile << durationValues[i] << "\n";
                }
                outputFile.close();
            } else {
                std::cerr << "Error opening output file.\n";
            }

    system("pause");
    return 0;
}
