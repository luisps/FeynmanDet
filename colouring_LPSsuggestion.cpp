//  Path_simulator
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include <string.h>
#include "layer.hpp"
#include "circuit.cpp"


// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
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
		//case 13:
		{
			out=2;
			break;
		}
	}
	return out;

}


// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
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
	//Erro=-1 
    
//map with colouring for each qubit at each intermediate state
    char RG_map_fw[NQ][L+1];
	char RG_map_bw[NQ][L+1];
	char RG_map[NQ][L+1];

// set all RG_map to RED
	for (int l = 0; l < L + 1; l++) {
		for (int nq = 0; nq < NQ; nq++) {
			RG_map_fw[nq][l] = 2;
			RG_map_bw[nq][l] = 2;
			RG_map[nq][l] = 2;
		}
	}


// set qubits in columns 0 and L according to initial and final_state
	for (int nq=0; nq<NQ;nq++){
	    RG_map_fw[nq][0] = init_state[nq];
	    RG_map_bw[nq][L] = final_state[nq];
	}	


// FS
// Iterate over the layers
for (int l = 0; l < L; l++) {
	for (int q = 0; q < NQ; q++) {
        RG_map_fw[q][l+1] = RG_map_fw[q][l];
    }
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
					//ver numero que esta no mapa para este qubit nesta layer: 0,1 ou 2
					input_RG=RG_map_fw[qubit][l];
					
					//se for vermelho, sai
					if (input_RG==2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					output_RG = gate_output_1(input_RG, name);
					
					// set colors in RG_map_fw
					RG_map_fw[qubit][l+1]=output_RG;
					
					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_fw[nq][l+1]==2){
    						continue;
						}
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
					//ver n�mero que est� no mapa para este qubit nesta layer: 0,1 ou 2
					input_RG=RG_map_fw[qubit][l];
					
					//se for vermelho, sai
					if (input_RG==2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					output_RG = gate_output_1(input_RG, name);
					
					// set colors in RG_map_fw
					RG_map_fw[qubit][l+1]=output_RG;


					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_fw[nq][l+1]==2){
    						continue;;
						}
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
					int input_RG_c;
					int input_RG_t;
					input_RG_c=RG_map_fw[c_qubit][l];
					input_RG_t=RG_map_fw[t_qubit][l];
                    input_RG[0] = input_RG_c;
                    input_RG[1] = input_RG_t;
					
				
					//se for vermelho, sai
					if (input_RG_c==2 && input_RG_t==2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					int output[2];
					gate_output_2(input_RG, name, output);
					// set colors in RG_map_fw
					RG_map_fw[c_qubit][l+1] = output[0];
					// Para cada qubit
					if (RG_map_fw[t_qubit][l] == 2) {
						// O qubit está num estado de superposição (2), então, não altere para 0 ou 1 sem verificar
						RG_map_fw[t_qubit][l+1] = 2;
					} else {
						// Normalmente aplica a lógica da CNOT
						RG_map_fw[t_qubit][l+1] = output[1];
					}

					
					//ver se o target e control ficam os 2 a vermelho
					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_fw[nq][l+1]==2){
    						continue;
						}
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
					int input_RG_c;
					int input_RG_t;
					input_RG_c=RG_map_fw[c_qubit][l];
					input_RG_t=RG_map_fw[t_qubit][l];
                    input_RG[0] = input_RG_c;
                    input_RG[1] = input_RG_t;
					
					
					
					//se for vermelho, sai
					if (input_RG_c==2 && input_RG_t==2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					int output[2];
					gate_output_2(input_RG, name, output);
					// set colors in RG_map_fw
					RG_map_fw[c_qubit][l+1] = output[0];
					// Para cada qubit
					if (RG_map_fw[t_qubit][l] == 2) {
						// O qubit está num estado de superposição (2), então, não altere para 0 ou 1 sem verificar
						RG_map_fw[t_qubit][l+1] = 2;
					} else {
						// Normalmente aplica a lógica da CNOT
						RG_map_fw[t_qubit][l+1] = output[1];
					}

					
					//ver se o target e control ficam os 2 a vermelho
					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_fw[nq][l+1]==2){
    						continue;
						}
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
/*
printf("FORWARD COLOURING (RG_map_fw):\n");
	for (int nq = 0; nq < NQ; nq++) {
		printf("Qubit %d: ", nq);
		for (int l = 0; l <= L; l++) {
			printf("%d ", RG_map_fw[nq][l]);
		}
		printf("\n");
	}
*/

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
					//ver n�mero que est� no mapa para este qubit nesta layer: 0,1 ou 2
					output_RG=RG_map_bw[qubit][l+1];
					
					//se for vermelho, sai
					if (output_RG==2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					input_RG = gate_output_1(output_RG, name);
					
					// set colors in RG_map_bw
					RG_map_bw[qubit][l]=input_RG;
					
					
					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_bw[nq][l+1]==2){
    						continue;
						}
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
					//ver n�mero que est� no mapa para este qubit nesta layer: 0,1 ou 2
					output_RG=RG_map_bw[qubit][l+1];
					
					//se for vermelho, sai
					if (output_RG==2) {
                        continue;
                    }
                    
					//compute output colour for this gate
					input_RG = gate_output_1(output_RG, name);
					
					// set colors in RG_map_bw
					RG_map_bw[qubit][l]=input_RG;
					
					
					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_bw[nq][l+1]==2){
    						continue;
						}
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

					// Cores de saída nesta layer (isto é, cores no tempo seguinte)
					int output_RG_c = RG_map_bw[c_qubit][l + 1];
					int output_RG_t = RG_map_bw[t_qubit][l + 1];

					// Se ambos forem 2 (irrelevantes), continua
					if (output_RG_c == 2 && output_RG_t == 2) {
						continue;
					}

					int output_RG[2] = {output_RG_c, output_RG_t};
					int input[2];
					gate_output_2(output_RG, name, input);  // Calcula input RG a partir de output RG

					// Atualiza RG_map_bw (estamos no backward sweep!)
					RG_map_bw[c_qubit][l] = input[0];

					if (output_RG_t == 2) {
						RG_map_bw[t_qubit][l] = 2;
					} else {
						RG_map_bw[t_qubit][l] = input[1];
					}

					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_fw[nq][l+1]==2){
    						continue;;
						}
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
			
					// Cores de saída nesta layer (estado depois da porta)
					int output_RG_c = RG_map_bw[c_qubit][l + 1];
					int output_RG_t = RG_map_bw[t_qubit][l + 1];
			
					// Se ambos forem irrelevantes (2), não há nada a fazer
					if (output_RG_c == 2 && output_RG_t == 2) {
						continue;
					}
			
					int output_RG[2] = {output_RG_c, output_RG_t};
					int input[2];
					gate_output_2(output_RG, name, input);  // Calcula input RG a partir de output RG
			
					// Atualiza RG_map_bw (backward sweep!)
					RG_map_bw[c_qubit][l] = input[0];
			
					if (output_RG_t == 2) {
						RG_map_bw[t_qubit][l] = 2;
					} else {
						RG_map_bw[t_qubit][l] = input[1];
					}

					for (int nq=0; nq<NQ;nq++){
    					if(RG_map_fw[nq][l+1]==2){
    						continue;;
						}
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
/*
printf("BACKWARD COLOURING (RG_map_bw):\n");
	for (int nq = 0; nq < NQ; nq++) {
		printf("Qubit %d: ", nq);
		for (int l = 0; l <= L; l++) {
			printf("%d ", RG_map_bw[nq][l]);
		}
		printf("\n");
	}           
*/          
	//ver a cor e binary digit final
	
	for (int a=0; a<=L; a++){
		for (int b=0; b<NQ; b++){
			if (RG_map_fw[b][a]==2 && RG_map_bw[b][a]==2){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=2;
			}
			
			else if (RG_map_fw[b][a]==2){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=RG_map_bw[b][a];
			}
			else if (RG_map_bw[b][a]==2){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=RG_map_fw[b][a];
			}
			
			else if (RG_map_fw[b][a]==RG_map_bw[b][a]){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=RG_map_fw[b][a];
			} else {
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=-1; //caso de conflito
			}
		}
		//printf("\n");
	
	}


	/*
	printf("COLOURING FINAL:\n");
	for (int nq = 0; nq < NQ; nq++) {
		printf("Qubit %d: ", nq);
		for (int l = 0; l <= L; l++) {
			printf("%d ", RG_map[nq][l]);
		}
		printf("\n");
	} 
*/
	int* final_states_output = new int[(L-1)*NQ];
	//guardar no formato de um unico array

	for (int q = 0; q < NQ; q++) {
		for (int l = 1; l < L; l++) {
			final_states_output[q * (L - 1) + (l - 1)] = RG_map[q][l];
		}
	}
	
	/*
	printf("FINAL STATES OUTPUT:\n");
for (int q = 0; q < NQ; q++) {
    printf("Qubit %d: ", q);
    for (int l = 0; l < L - 1; l++) {
        printf("%d ", final_states_output[q * (L - 1) + l]);
    }
    printf("\n");
}
*/
    return final_states_output;
}