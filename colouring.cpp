//  Path_simulator
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include <string.h>
#include "layer.hpp"
#include "circuit.cpp"


// pega no input estado inicial e gate -> da o estado final (para ver se e 1 ou 0)
//gates 1 qubit
int gate_output_1(int init_state, int name){
	int out=0;	
	switch (name){
		case 0:
		case 4:
		case 5:
		case 6:
        case 13:
		{ //id
			if(init_state==1){
				out=1;
			}
			break;
		}
		case 2:
		case 3:
		{ //X //Y
			if(init_state==0){
				out=1;
			}
			break;
		}
		
		case 1:
		case 11:
		case 12:
		{ //RX
			out=2;
			break;
		}
	}
	return out;
}


// pega no input estado inicial e gate -> da o estado final (para ver se é 1 ou 0)
// gates 2 qubits
void gate_output_2(int* init_state, int name, int* out){
	out[0]=0;
	out[1]=0;
	switch (name){
		case 20:
		case 22:
		case 31: { 
			if(init_state[0]==1){
				out[0]=1;
			}
			if(init_state[1]==1){
				out[1]=1;
			}
			break;
		}
		case 21: { //CX
			if(init_state[0]==1){
				out[0]=1;
				if(init_state[1]==0){
					out[1]=1;
				}	
			} else{ 
				if(init_state[1]==1){
					out[1]=1;
				}
			}
			break;
		}
	}
}


int* fs_bs_RG(TCircuit *circuit, int* init_state, int* final_state) {
    // Layers and qubits
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits;
    const int N = (int)(powf(2.f, (float)NQ));

    // Colors
    // RED==2;
    // GREEN_N0==0;
    // GREEN_N1==1;
    // Error=-1 

    // Map with colouring for each qubit at each intermediate state
    char RG_map_fw[NQ][L + 1];
    char RG_map_bw[NQ][L + 1];
    char RG_map[NQ][L + 1];

    // Set all RG_map to RED (2)
    for (int l = 0; l < L + 1; l++) {
        for (int nq = 0; nq < NQ; nq++) {
            RG_map_fw[nq][l] = 2;
            RG_map_bw[nq][l] = 2;
        }
    }

    // Set qubits in columns 0 and L according to initial and final state
    for (int nq = 0; nq < NQ; nq++) {
        RG_map_fw[nq][0] = init_state[nq];
        RG_map_bw[nq][L] = final_state[nq];
    }

    // Forward simulation (FS)
    for (int l = 0; l < L; l++) {
        TCircuitLayer* layer = &circuit->layers[l];

        for (int GT = 0; GT < 4; GT++) {
            void* gates_ptr = layer->gates[GT];
            int num_gates = layer->num_type_gates[GT];

            switch (GT) {
                case 0: {
                    TGate1P0* g1p0_ptr = (TGate1P0*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g1p0_ptr++) {
                        int name = g1p0_ptr->name;
                        int qubit = g1p0_ptr->qubit;
                        int input_RG = RG_map_fw[qubit][l];
                        int output_RG;

                        if (input_RG == 2) continue;  // Skip if RED (2)

                        output_RG = gate_output_1(input_RG, name);
                        RG_map_fw[qubit][l + 1] = output_RG;
                    }
                    break;
                }

                case 1: {
                    TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                        int name = g1p1_ptr->fdata.name;
                        int qubit = g1p1_ptr->fdata.qubit;
                        int input_RG = RG_map_fw[qubit][l];
                        int output_RG;

                        if (input_RG == 2) continue;  // Skip if RED (2)

                        output_RG = gate_output_1(input_RG, name);
                        RG_map_fw[qubit][l + 1] = output_RG;
                    }
                    break;
                }

                case 2: {
                    TGate2P0* g2p0_ptr = (TGate2P0*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g2p0_ptr++) {
                        int name = g2p0_ptr->name;
                        int c_qubit = g2p0_ptr->c_qubit;
                        int t_qubit = g2p0_ptr->t_qubit;

                        int input_RG[2];
                        input_RG[0] = RG_map_fw[c_qubit][l];
                        input_RG[1] = RG_map_fw[t_qubit][l];

                        if (input_RG[0] == 2 && input_RG[1] == 2) continue;  // Skip if both RED (2)

                        int output[2];
                        gate_output_2(input_RG, name, output);

                        RG_map_fw[t_qubit][l + 1] = output[1];
                    }
                    break;
                }

                case 3: {
                    TGate2P1* g2p1_ptr = (TGate2P1*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g2p1_ptr++) {
                        int name = g2p1_ptr->fdata.name;
                        int c_qubit = g2p1_ptr->fdata.c_qubit;
                        int t_qubit = g2p1_ptr->fdata.t_qubit;

                        int input_RG[2];
                        input_RG[0] = RG_map_fw[c_qubit][l];
                        input_RG[1] = RG_map_fw[t_qubit][l];

                        if (input_RG[0] == 2 && input_RG[1] == 2) continue;  // Skip if both RED (2)

                        int output[2];
                        gate_output_2(input_RG, name, output);

                        RG_map_fw[t_qubit][l + 1] = output[1];
                    }
                    break;
                }

                default:
                    printf("Gate type invalid");
                    break;
            }
        }
    }

    // Backward simulation (BS)
    for (int l = L - 1; l >= 0; l--) {
        TCircuitLayer* layer = &circuit->layers[l];

        for (int GT = 0; GT < 4; GT++) {
            void* gates_ptr = layer->gates[GT];
            int num_gates = layer->num_type_gates[GT];

            switch (GT) {
                case 0: {
                    TGate1P0* g1p0_ptr = (TGate1P0*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g1p0_ptr++) {
                        int name = g1p0_ptr->name;
                        int qubit = g1p0_ptr->qubit;
                        int output_RG = RG_map_bw[qubit][l + 1];
                        int input_RG;

                        if (output_RG == 2) continue;  // Skip if RED (2)

                        input_RG = gate_output_1(output_RG, name);
                        RG_map_bw[qubit][l] = input_RG;
                    }
                    break;
                }

                case 1: {
                    TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                        int name = g1p1_ptr->fdata.name;
                        int qubit = g1p1_ptr->fdata.qubit;
                        int output_RG = RG_map_bw[qubit][l + 1];
                        int input_RG;

                        if (output_RG == 2) continue;  // Skip if RED (2)

                        input_RG = gate_output_1(output_RG, name);
                        RG_map_bw[qubit][l] = input_RG;
                    }
                    break;
                }

                case 2: {
                    TGate2P0* g2p0_ptr = (TGate2P0*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g2p0_ptr++) {
                        int name = g2p0_ptr->name;
                        int c_qubit = g2p0_ptr->c_qubit;
                        int t_qubit = g2p0_ptr->t_qubit;

                        int output_RG[2];
                        output_RG[0] = RG_map_bw[c_qubit][l];
                        output_RG[1] = RG_map_bw[t_qubit][l];

                        if (output_RG[0] == 2 && output_RG[1] == 2) continue;  // Skip if both RED (2)

                        int input[2];
                        gate_output_2(output_RG, name, input);

                        RG_map_bw[t_qubit][l] = input[1];
                    }
                    break;
                }

                case 3: {
                    TGate2P1* g2p1_ptr = (TGate2P1*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g2p1_ptr++) {
                        int name = g2p1_ptr->fdata.name;
                        int c_qubit = g2p1_ptr->fdata.c_qubit;
                        int t_qubit = g2p1_ptr->fdata.t_qubit;

                        int output_RG[2];
                        output_RG[0] = RG_map_bw[c_qubit][l];
                        output_RG[1] = RG_map_bw[t_qubit][l];

                        if (output_RG[0] == 2 && output_RG[1] == 2) continue;  // Skip if both RED (2)

                        int input[2];
                        gate_output_2(output_RG, name, input);

                        RG_map_bw[t_qubit][l] = input[1];
                    }
                    break;
                }

                default:
                    printf("Gate type invalid");
                    break;
            }
        }
    }

   //ver a cor e binary digit final
	for (int b=0; b<NQ; b++){
		for (int a=0; a<L-1; a++){
			if (RG_map_fw[b][a]==2 && RG_map_bw[b][a]==2){
				RG_map[b][a]=2;
			}
			
			else if (RG_map_fw[b][a]==2){
				RG_map[b][a]=RG_map_bw[b][a];
			}
			else if (RG_map_bw[b][a]==2){
				RG_map[b][a]=RG_map_fw[b][a];
			}
			
			else if (RG_map_fw[b][a]==RG_map_bw[b][a]){
					RG_map[b][a]=RG_map_fw[b][a];
			} else {
				RG_map[b][a]=-1; //caso de conflito
			}
		}
		//printf("\n");
	
	}

	int* final_states_output = new int[(L-1)*NQ];
	//guardar no formato de um unico array
	
	for (int q=0; q<NQ;q++){
    	for (int l=0; l<L-1; l++){
    		final_states_output[q*(L-1)+l]=RG_map[q][l];
		}
	}
    return final_states_output;
}


