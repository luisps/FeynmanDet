//
//  simulation_attempt.cpp
//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//


#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include <string.h>
#include "layer.hpp"
#include "circuit.cpp"


// pega no input estado inicial e gate -> dá o estado final (para ver se é 1 ou 0)
//gates 1 qubit
int gate_output_1(int init_state, int name){
	int out=0;
	
	switch (name){
		case 0:
		case 4:
		case 5:
		case 6:
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
		case 13:
		{ //RX
			out=2;
			break;
		}
	}
	return out;

}


// pega no input estado inicial e gate -> dá o estado final (para ver se é 1 ou 0)
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



int* fs_bs_RG(TCircuit *circuit, int* init_state,  int* final_state){	
	
//Layers e qubits    
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits; 
    const int N = (int)(powf(2.f, (float)NQ));
    
//Colours
	//RED==2;
	//GREEN_N0==0;
	//GREEN_N1==1;    
    
//map with colouring for each qubit at each intermediate state

    char RG_map_fw[NQ][L+1];
	char RG_map_bw[NQ][L+1];
	char RG_map[NQ][L+1];

// set all RG_map to RED

//Using for loops:
/*    for (int l=0; l<L+1; l++){
		for (int nq=0; nq<NQ;nq++){
    		RG_map_fw[l+1][nq] = 2;
    		RG_map_bw[l+1][nq] = 2;
		}	
	}
*/

//Using memset:

// Set all elements of to 2 (red)
	memset(RG_map_fw, 2, (L + 1) * NQ * sizeof(char));
	memset(RG_map_bw, 2, (L + 1) * NQ * sizeof(char));


// set qubits in columns 0 and L to the corresponding 
//GREEN0 or GREEN1 according to initial_state and final_state
	
//RG_map_fw[…][0] set <- according to initial state
//RG_map_bw[…][L] set <- according to final state

	for (int nq=0; nq<NQ;nq++){
	    RG_map_fw[nq][0] = initial_state[nq];
	    RG_map_bw[nq][L] = final_state[nq];
	}	
		

// FS

// Iterate over the layers
for (int l = 0; l < L; l++) {

//layer <- get_layer_from_circuit (circuit, l)
    TCircuitLayer* layer = &circuit->layers[l]; 

// Iterate over the gate types for each layer
    for (int GT = 0; GT < 4; GT++) { 
        void* gates_ptr = layer->gates[GT];
        int num_gates = layer->num_type_gates[GT];
        
        switch (GT) {
        	
            case 0: {
                TGate1P0* g1p0_ptr = (TGate1P0*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p0_ptr++) {
                    int name = g1p0_ptr->name;
                    int qubit = g1p0_ptr->qubit; //qubit que participa na gate
					int input_RG;
					int output_RG;
					
					//get colors from RG_map_fw
					//ver número que está no mapa para este qubit nesta layer: 0,1 ou 2
					input_RG=RG_map_fw[qubit][l];
					
					//se for vermelho, sai
					if (input_RG=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					output_RG = gate_output_1(input_RG, name);
					
					// set colors in RG_map_fw
					RG_map_fw[qubit][l+1]=output_RG;
					
					
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_fw[nq][l+1] = 2;
						}){ break;
						}
				}
				break;
			}
				
				
			case 1: {
                TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                    int name = g1p1_ptr->fdata.name;
                    int qubit = g1p1_ptr->fdata.qubit;
                    int input_RG;
					int output_RG;
					//get colors from RG_map_fw
					//ver número que está no mapa para este qubit nesta layer: 0,1 ou 2
					input_RG=RG_map_fw[qubit][l];
					
					//se for vermelho, sai
					if (input_RG=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					output_RG = gate_output_1(input_RG, name);
					
					// set colors in RG_map_fw
					RG_map_fw[qubit][l+1]=output_RG;
					
					
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_fw[nq][l+1] = 2;
						}){ break;
						}
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
                    input_RG[0] = input_RG_c;
                    input_RG[1] = input_RG_t;
					
					input_RG_c=RG_map_fw[c_qubit][l];
					input_RG_t=RG_map_fw[t_qubit][l];
					
					//se for vermelho, sai
					if (input_RG_c=2 && input_RG_t=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					int output[2];
					gate_output_2(input_RG, name, output);
					// set colors in RG_map_fw
					RG_map_fw[t_qubit][l+1]=output[1];
					
					//ver se o target e control ficam os 2 a vermelho
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_fw[nq][l] = 2;
						}){ break;
						}
					
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
                    input_RG[0] = input_RG_c;
                    input_RG[1] = input_RG_t;
					
					input_RG_c=RG_map_fw[c_qubit][l];
					input_RG_t=RG_map_fw[t_qubit][l];
					
					//se for vermelho, sai
					if (input_RG_c=2 && input_RG_t=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					int output[2];
					gate_output_2(input_RG, name, output);
					// set colors in RG_map_fw
					RG_map_fw[t_qubit][l+1]=output[1];
					
					//ver se o target e control ficam os 2 a vermelho
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_fw[nq][l] = 2;
						}){ break;
						}
					
				}
                break;
            }		
			default: {
				printf("Gate type invalid");
                break;	
			}
		}
	}
}




// BS
for (int l = L - 1; l >= 0; l--){ // Iterate over the layers
    TCircuitLayer* layer = &circuit->layers[l];

    for (int GT = 0; GT < 4; GT++){ // Iterate over the gate types
        void* gates_ptr = layer->gates[GT];
        int num_gates = layer->num_type_gates[GT];
        
        switch (GT) {
            case 0: {
                TGate1P0* g1p0_ptr = (TGate1P0*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p0_ptr++) {
                    int name = g1p0_ptr->name;
                    int qubit = g1p0_ptr->qubit;
				    int input_RG;
					int output_RG;                    
					
					//get colors from RG_map_bw
					//ver número que está no mapa para este qubit nesta layer: 0,1 ou 2
					output_RG=RG_map_bw[qubit][l+1];
					
					//se for vermelho, sai
					if (output_RG=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					input_RG = gate_output_1(output_RG, name);
					
					// set colors in RG_map_fw
					RG_map_bw[qubit][l]=input_RG;
					
					
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_bw[nq][l] = 2;
						}){ break;
						}
				}
				break;
			}
			
            case 1: {
                TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                    int name = g1p1_ptr->fdata.name;
                    int qubit = g1p1_ptr->fdata.qubit;
					int input_RG;
					int output_RG;
					//get colors from RG_map_bw
					//ver número que está no mapa para este qubit nesta layer: 0,1 ou 2
					output_RG=RG_map_bw[qubit][l+1];
					
					//se for vermelho, sai
					if (output_RG=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					input_RG = gate_output_1(output_RG, name);
					
					// set colors in RG_map_fw
					RG_map_bw[qubit][l]=input_RG;
					
					
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_bw[nq][l] = 2;
						}){ break;
						}
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
                    output_RG[0] = input_RG_c;
                    output_RG[1] = input_RG_t;
					
					output_RG_c=RG_map_bw[c_qubit][l];
					output_RG_t=RG_map_bw[t_qubit][l];
					
					//se for vermelho, sai
					if (output_RG_c=2 && output_RG_t=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					int input[2];
					gate_output_2(output_RG, name, input);
					// set colors in RG_map_fw
					RG_map_fw[t_qubit][l]=input[1];
					
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_bw[nq][l] = 2;
						}){ break;
						}
					
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
                    output_RG[0] = input_RG_c;
                    output_RG[1] = input_RG_t;
					
					output_RG_c=RG_map_bw[c_qubit][l];
					output_RG_t=RG_map_bw[t_qubit][l];
					
					//se for vermelho, sai
					if (output_RG_c=2 && output_RG_t=2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					int input[2];
					gate_output_2(output_RG, name, input);
					// set colors in RG_map_fw
					RG_map_fw[t_qubit][l]=input[1];
					
					if(for (int nq=0; nq<NQ;nq++){
    					RG_map_bw[nq][l] = 2;
						}){ break;
						}
					
				}
                break;
            }  					

            default: {
				printf("Gate type invalid");
                break;
            
			}
        }
    }
}            
            
	printf("\n This is the forward sweep \n");
    for (int q=0; q<NQ;q++){
    	printf("%d ", pos_f[q]);
		printf("\n");
	}
	printf("\n This is backward sweep \n");
    for (int q=0; q<NQ;q++){
    	printf("%d ", pos[q]);
		printf("\n");
	}

	
	printf("This is the full sweep \n");
	
//?????
	
	
    return final_states_output;
}









