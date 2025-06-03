//
//  simulation_attempt.cpp
//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//


#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
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
	const int green=0;
	const int red=1;
 //nomear as branching gates de acordo com tipo GP
    int namesBG1P0[1] = {1}; //H
    int namesBG1P1[3] = {11, 12, 13}; //rx,ry,rz
    int namesBG2P0[1] = {21}; //CX
    int *namesB[3] = {namesBG1P0,namesBG1P1,namesBG2P0};
    
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits; 
    const int N = (int)(powf(2.f, (float)NQ));

 
int first_branching_gate_encountered[NQ];
for (int q=0; q<NQ;q++){
    	first_branching_gate_encountered[q] = 0;
}
int pos[NQ];
int inter_states[L-1][NQ];
int previous_state;

// BS first only
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
                    int branching = (name == namesB[0][0]) ? 1 : 0;
                    int qubit = g1p0_ptr->qubit;
					
					//COLOR
                    // Check if this is the first branching gate for the qubit
                    if (branching == 1 && first_branching_gate_encountered[qubit]!=1){
                        first_branching_gate_encountered[qubit] = 1;
                        pos[qubit] = l-1;
                    }
                    
                    //Check if it is the last gate
					if (l == 0 && first_branching_gate_encountered[qubit]!=1){
                    	pos[qubit] = -1;
					}
					
					//STATE
                    //If it is the first gate, the previous state is the initial state
                    if (l == L-1) {
                        previous_state = final_state[qubit];
                    } else {
                        previous_state = inter_states[l-1][qubit];
                    }
                    
                    //Compute state from the previous
                    if (first_branching_gate_encountered[qubit]!=1) {
                        inter_states[l][qubit] = gate_output_1(previous_state, name);
                    } else {
                        inter_states[l][qubit] = 2;
                    }
                }
                break;
            }
            case 1: {
                TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                    int name = g1p1_ptr->fdata.name;
                    bool branching = (name == namesB[1][0]) ? 1 : 0;
                    branching += ((name == namesB[1][1]) ? 1 : 0);
                    branching += ((name == namesB[1][2]) ? 1 : 0);
                    int qubit = g1p1_ptr->fdata.qubit;

                    
					//COLOR
					if (branching == 1 && first_branching_gate_encountered[qubit]!=1) {
                        first_branching_gate_encountered[qubit] = 1;
                        pos[qubit] = l-1;   
                    }
					
					if (l == 0 && first_branching_gate_encountered[qubit]!=1) {
                    	pos[qubit] = -1;
					}
                    
					//STATE
					if (l == L-1) {
                        previous_state = final_state[qubit];
                    } else {
                        previous_state = inter_states[l-1][qubit];
                    }
                    
                    if (first_branching_gate_encountered[qubit]!=1) {
                        inter_states[l][qubit] = gate_output_1(previous_state, name);
                    } else {
                        inter_states[l][qubit] = 2;
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
					
					
					//COLOR
                    if (first_branching_gate_encountered[c_qubit]==1 && first_branching_gate_encountered[t_qubit]!=1) {
                        first_branching_gate_encountered[t_qubit] = 1;
                        pos[t_qubit] = l-1;
                        
                    } 
					
					if (l == 0 && first_branching_gate_encountered[t_qubit]!=1) {
                    	pos[t_qubit] = -1;
                    	
					}
					if (l == 0 && first_branching_gate_encountered[c_qubit]!=1) {
                    	pos[c_qubit] = -1;
                    
					}
                    
                    
                    //STATE
                    int previous_state_t;
                    int previous_state_c;
                    
                    if (l == L-1) {
                        previous_state_c = final_state[c_qubit];
                        previous_state_t = final_state[t_qubit];
                    } else {
                        previous_state_c = inter_states[l-1][c_qubit];
                        previous_state_t = inter_states[l-1][t_qubit];
                    }
                    
                    int previous_aux[2];
                    previous_aux[0] = previous_state_c;
                    previous_aux[1] = previous_state_t;
                    
                    if (first_branching_gate_encountered[c_qubit]==1) {
                        inter_states[l][c_qubit] = 2;
                        inter_states[l][t_qubit] = 2;
                    } else if (first_branching_gate_encountered[t_qubit]==1 && first_branching_gate_encountered[c_qubit]!=1) {
                        inter_states[l][t_qubit] = 2;
                        inter_states[l][c_qubit] = inter_states[l-1][c_qubit];
                    } else {
                        int output[2];
                        gate_output_2(previous_aux, name, output);
                        inter_states[l][t_qubit] = output[1];
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

					//COLOR
					if (first_branching_gate_encountered[c_qubit]==1 && first_branching_gate_encountered[t_qubit]!=1) {
                        first_branching_gate_encountered[t_qubit] = 1;
                        pos[t_qubit] = l-1;
                       
                    } 
                    
					if (l == 0 && first_branching_gate_encountered[t_qubit]!=1) {
                    	pos[t_qubit] = -1;
                    
					}
					
					if (l == 0 && !first_branching_gate_encountered[c_qubit]) {
                    	pos[c_qubit] = -1;
                    
                	}
                    
                    
                    //STATE
                    int previous_state_t;
                    int previous_state_c;
                    
                    if (l == L-1) {
                        previous_state_c = final_state[c_qubit];
                        previous_state_t = final_state[t_qubit];
                    } else {
                        previous_state_c = inter_states[l-1][c_qubit];
                        previous_state_t = inter_states[l-1][t_qubit];
                    }
                    
                    int previous_aux[2];
                    previous_aux[0] = previous_state_c;
                    previous_aux[1] = previous_state_t;
                    
                    if (first_branching_gate_encountered[c_qubit]==1) {
                        inter_states[l][c_qubit] = 2;
                        inter_states[l][t_qubit] = 2;
                    } else if (first_branching_gate_encountered[t_qubit]==1 && first_branching_gate_encountered[c_qubit]!=1) {
                        inter_states[l][t_qubit] = 2;
                        inter_states[l][c_qubit] = inter_states[l-1][c_qubit];
                    } else {
                        int output[2];
                        gate_output_2(previous_aux, name, output);
                        inter_states[l][t_qubit] = output[1];
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




int first_branching_gate_encountered_f[NQ];
int inter_states_f[L-1][NQ];
int pos_f[NQ];

for (int q=0; q<NQ;q++){
    	first_branching_gate_encountered_f[q] = 0;
}


// FS first only
for (int l = 0; l < L; l++) { // Iterate over the layers
    TCircuitLayer* layer = &circuit->layers[l];
    for (int GT = 0; GT < 4; GT++) { // Iterate over the gate types
        void* gates_ptr = layer->gates[GT];
        int num_gates = layer->num_type_gates[GT];
        
        switch (GT) {
            case 0: {
                TGate1P0* g1p0_ptr = (TGate1P0*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p0_ptr++) {
                    int name = g1p0_ptr->name;
                    int branching = (name == namesB[0][0]) ? 1 : 0;
                    int qubit = g1p0_ptr->qubit;
					
					//COLOR
                    // Check if this is the first branching gate for the qubit
                    if (branching == 1 && first_branching_gate_encountered_f[qubit]!=1){
                        first_branching_gate_encountered_f[qubit] = 1;
                        pos_f[qubit] = l;
                    }
                    
                    //Check if it is the last gate
					if (l == L-1 && first_branching_gate_encountered_f[qubit]!=1){
                    	pos_f[qubit] = L-1;
                  
					}
					
					//STATE
                    //If it is the first gate, the previous state is the initial state
                    if (l == 0) {
                        previous_state = init_state[qubit];
                    } else {
                        previous_state = inter_states_f[l-1][qubit];
                    }
                    
                    //Compute state from the previous
                    if (first_branching_gate_encountered_f[qubit]!=1) {
                        inter_states_f[l][qubit] = gate_output_1(previous_state, name);
                    } else {
                        inter_states_f[l][qubit] = 2;
                    }
                }
                break;
            }
            case 1: {
                TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                    int name = g1p1_ptr->fdata.name;
                    bool branching = (name == namesB[1][0]) ? 1 : 0;
                    branching += ((name == namesB[1][1]) ? 1 : 0);
                    branching += ((name == namesB[1][2]) ? 1 : 0);
                    int qubit = g1p1_ptr->fdata.qubit;

                    
					//COLOR
					if (branching == 1 && first_branching_gate_encountered_f[qubit]!=1) {
                        first_branching_gate_encountered_f[qubit] = 1;
                        pos_f[qubit] = l;
                        
                    }
					
					if (l == L-1 && first_branching_gate_encountered_f[qubit]!=1) {
                    	pos_f[qubit] = L-1;
                    
					}
                    
					//STATE
					if (l == 0) {
                        previous_state = init_state[qubit];
                    } else {
                        previous_state = inter_states_f[l-1][qubit];
                    }
                    
                    if (first_branching_gate_encountered_f[qubit]!=1) {
                        inter_states_f[l][qubit] = gate_output_1(previous_state, name);
                    } else {
                        inter_states_f[l][qubit] = 2;
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
					
					
					//COLOR
                    if (first_branching_gate_encountered_f[c_qubit]==1 && first_branching_gate_encountered_f[t_qubit]!=1) {
                        first_branching_gate_encountered_f[t_qubit] = 1;
                        pos_f[t_qubit] = l;
                        
                    } 
					
					if (l == L-1 && first_branching_gate_encountered_f[t_qubit]!=1) {
                    	pos_f[t_qubit] = L-1;
                    	
					}
					if (l == L-1 && first_branching_gate_encountered_f[c_qubit]!=1) {
                    	pos_f[c_qubit] = L-1;
                    
					}
                    
                    
                    //STATE
                    int previous_state_t;
                    int previous_state_c;
                    
                    if (l == 0) {
                        previous_state_c = init_state[c_qubit];
                        previous_state_t = init_state[t_qubit];
                    } else {
                        previous_state_c = inter_states_f[l-1][c_qubit];
                        previous_state_t = inter_states_f[l-1][t_qubit];
                    }
                    
                    int previous_aux[2];
                    previous_aux[0] = previous_state_c;
                    previous_aux[1] = previous_state_t;
                    
                    if (first_branching_gate_encountered_f[c_qubit]==1) {
                        inter_states_f[l][c_qubit] = 2;
                        inter_states_f[l][t_qubit] = 2;
                    } else if (first_branching_gate_encountered_f[t_qubit]==1 && first_branching_gate_encountered_f[c_qubit]!=1) {
                        inter_states_f[l][t_qubit] = 2;
                        inter_states_f[l][c_qubit] = inter_states_f[l-1][c_qubit];
                    } else {
                        int output[2];
                        gate_output_2(previous_aux, name, output);
                        inter_states_f[l][t_qubit] = output[1];
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

					//COLOR
					if (first_branching_gate_encountered_f[c_qubit]==1 && first_branching_gate_encountered_f[t_qubit]!=1) {
                        first_branching_gate_encountered_f[t_qubit] = 1;
                        pos_f[t_qubit] = l;
                       
                    } 
                    
					if (l == L-1 && first_branching_gate_encountered_f[t_qubit]!=1) {
                    	pos_f[t_qubit] = L-1;
                    
					}
					
					if (l == L-1 && !first_branching_gate_encountered_f[c_qubit]) {
                    	pos_f[c_qubit] = L-1;
                	}
                    
                    
                    //STATE
                    int previous_state_t;
                    int previous_state_c;
                    
                    if (l == 0) {
                        previous_state_c = init_state[c_qubit];
                        previous_state_t = init_state[t_qubit];
                    } else {
                        previous_state_c = inter_states_f[l-1][c_qubit];
                        previous_state_t = inter_states_f[l-1][t_qubit];
                    }
                    
                    int previous_aux[2];
                    previous_aux[0] = previous_state_c;
                    previous_aux[1] = previous_state_t;
                    
                    if (first_branching_gate_encountered_f[c_qubit]==1) {
                        inter_states_f[l][c_qubit] = 2;
                        inter_states_f[l][t_qubit] = 2;
                    } else if (first_branching_gate_encountered_f[t_qubit]==1 && first_branching_gate_encountered_f[c_qubit]!=1) {
                        inter_states_f[l][t_qubit] = 2;
                        inter_states_f[l][c_qubit] = inter_states_f[l-1][c_qubit];
                    } else {
                        int output[2];
                        gate_output_2(previous_aux, name, output);
                        inter_states_f[l][t_qubit] = output[1];
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

	
	int* final_colour_output = new int[(L-1)*NQ];
	int final_states[L-1][NQ];
	int final_colour[L-1][NQ];
	
	printf("This is the full sweep \n");
	
	//ver a cor final (um vez o outro, pois 1x0 dá 0: ou seja, se for green num fica green)
	for (int b=0; b<NQ; b++){
		for (int a=0; a<L-1; a++){
			final_colour_output[a*(L-1)+b] = (a <= pos[b]) * (a >= pos_f[b]);
			final_colour[a][b] = (a <= pos[b]) * (a >= pos_f[b]);
			printf("%d ", final_colour[a][b]);
		}
		printf("\n");
	
	}
	

	//atribuir aos estados verdes o seu valor (1 ou 0)
	// se for verdes mas diferentes é impossível
	
	for(int l=0; l<L-1;l++){
		for(int q=0; q<NQ; q++){
			if(final_colour[l][q]==green){
				final_states[l][q]=inter_states[l][q];
			}
			//if((l > pos[q]) && (l < pos_f[q]) && inter_states_f[l][q] != inter_states[l][q]){
			//	final_states[l][q]=-1;
			//}
			if(final_colour[l][q]==red){
				final_states[l][q]=2;
			}
			
		}
	}
	
	int* final_states_output = new int[(L-1)*NQ];
	//guardar no formato de um único array
	
	//printf("These are the final states \n");
	for (int q=0; q<NQ;q++){
    	for (int l=0; l<L-1; l++){
    		final_states_output[q*(L-1)+l]=final_states[l][q];
    		//printf("%d ", final_states[l][q]);
		}
		//printf("\n");
	}
	
	//printf("Successfuly calculated states and colours \n");
	
	
    return final_states_output;
}








