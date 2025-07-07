//  Path_simulator
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
#include "colouring_RG.hpp"

#include <math.h>
#include <string.h>
#include <stdio.h>

// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
//gates 1 qubit

// Correspondence between gate names numbers and strings
// G1P0 - 1 qubit, no parameters
// 'id' - 0
// 'h'  - 1 B
// 'x'  - 2
// 'y'  - 3
// 'z'  - 4
// 's'  - 5
// 't'  - 6
// G1P1 - 1 qubit, 1 parameter
// 'rx'  - 11 B
// 'ry'  - 12 B
// 'rz'  - 13
// 'p'  - 14
static colourT gate_output_1(colourT init_state, int name){
    colourT out=GREEN0;
	
    if (init_state==RED) {
        out=RED;
    }
    else {
        switch (name){
            case 0:
            case 4:
            case 5:
            case 6:
            case 13:
            { //id
                if(init_state==GREEN1){
                    out=GREEN1;
                }
                break;
            }
                
            case 2:
            case 3:
            { //X //Y
                if(init_state==GREEN0){
                    out=GREEN1;
                }
                break;
            }
                
            case 1:
            case 11:
            case 12:
                //case 13:
            {
                out=RED;
                break;
            }
        }
    }
	return out;

}


// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
// gates 2 qubits
// G2P0 - 2 qubits, no parameters
// 'id2'  - 20    Used for errors
// 'cx'  - 21 B
// 'cz'  - 22
// G2P1 - 2 qubits, 1 parameter
// 'cp'  - 31

static void gate_output_2(colourT* init_state, int name, colourT* out){
	
	out[0]=GREEN0;
	out[1]=GREEN0;
    switch (name){
        case 20:
        case 22:
        case 31: {
            out[0]=init_state[0];
            out[1]=init_state[1];
            break;
        }
            
        case 21: { //CX
            
            if (init_state[0]==RED) { // ctrl ==2, trgt don't care
                out[0]=out[1]=RED;
            }
            else if (init_state[1]==RED) { // trgt ==2, ctrl !=2
                out[0]=init_state[0];
                out[1]=RED;
            }
            else if(init_state[0]==GREEN1){
                out[0]=GREEN1;
                if(init_state[1]==GREEN0){
                    out[1]=GREEN1;
                }
            } else{
                if(init_state[1]==GREEN1){
                    out[1]=GREEN1;
                }
            }
            
            break;
        }
    }
}


colourT* fs_bs_RG(TCircuit *circuit, int* init_state_arr,  int* final_state_arr){
	
//Layers e qubits    
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits; 
    
//Colours
	//RED==2;
	//GREEN0==0;
	//GREEN1==1;
	//Erro=-1 
    
//map with colouring for each qubit at each intermediate state
    colourT RG_map_fw[NQ][L+1];
    colourT RG_map_bw[NQ][L+1];
    colourT RG_map[NQ][L+1];

// set all RG_map to RED
	for (int l = 0; l < L + 1; l++) {
		for (int nq = 0; nq < NQ; nq++) {
            RG_map_fw[nq][l] = RED;
            RG_map_bw[nq][l] = RED;
            RG_map[nq][l] = RED;
		}
	}


// set qubits in columns 0 and L according to initial and final_state
	for (int nq=0; nq<NQ;nq++){
        RG_map_fw[nq][0] = (colourT)init_state_arr[nq];
        RG_map_bw[nq][L] = (colourT)final_state_arr[nq];
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
                        colourT input_RG;
                        colourT output_RG;
                        
                        //get colors from RG_map_fw
                        //ver numero que esta no mapa para este qubit nesta layer: 0,1 ou 2
                        input_RG=RG_map_fw[qubit][l];
                        
                        //se for vermelho, sai
                        if (input_RG==RED) {
                            continue;
                        }
                        
                        //compute output colour for this gate
                        output_RG = gate_output_1(input_RG, name);
                        
                        // set colors in RG_map_fw
                        RG_map_fw[qubit][l+1]=output_RG;
                        
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_fw[nq][l+1]==RED){
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
                        colourT input_RG;
                        colourT output_RG;
                        //get colors from RG_map_fw
                        //ver n�mero que est� no mapa para este qubit nesta layer: 0,1 ou 2
                        input_RG=RG_map_fw[qubit][l];
                        
                        //se for vermelho, sai
                        if (input_RG==RED) {
                            continue;
                        }
                        
                        //compute output colour for this gate
                        output_RG = gate_output_1(input_RG, name);
                        
                        // set colors in RG_map_fw
                        RG_map_fw[qubit][l+1]=output_RG;
                        
                        
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_fw[nq][l+1]==RED){
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
                        
                        
                        colourT input_RG[2];
                        colourT input_RG_c;
                        colourT input_RG_t;
                        input_RG_c=RG_map_fw[c_qubit][l];
                        input_RG_t=RG_map_fw[t_qubit][l];
                        input_RG[0] = input_RG_c;
                        input_RG[1] = input_RG_t;
                        
                        
                        //se for vermelho, sai
                        if (input_RG_c==RED && input_RG_t==RED) {
                            continue;
                        }
                        
                        //compute output colour for this gate
                        colourT output[2];
                        gate_output_2(input_RG, name, output);
                        // set colors in RG_map_fw
                        RG_map_fw[c_qubit][l+1] = output[0];
                        // Para cada qubit
                        if (RG_map_fw[t_qubit][l] == RED) {
                            // O qubit está num estado de superposição (2), então, não altere para 0 ou 1 sem verificar
                            RG_map_fw[t_qubit][l+1] = RED;
                        } else {
                            // Normalmente aplica a lógica da CNOT
                            RG_map_fw[t_qubit][l+1] = output[1];
                        }
                        
                        
                        //ver se o target e control ficam os 2 a vermelho
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_fw[nq][l+1]==RED){
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
                        
                        colourT input_RG[2];
                        colourT input_RG_c;
                        colourT input_RG_t;
                        input_RG_c=RG_map_fw[c_qubit][l];
                        input_RG_t=RG_map_fw[t_qubit][l];
                        input_RG[0] = input_RG_c;
                        input_RG[1] = input_RG_t;
                        
                        //se for vermelho, sai
                        if (input_RG_c==RED && input_RG_t==RED) {
                            continue;
                        }
                        
                        //compute output colour for this gate
                        colourT output[2];
                        gate_output_2(input_RG, name, output);
                        // set colors in RG_map_fw
                        RG_map_fw[c_qubit][l+1] = output[0];
                        // Para cada qubit
                        if (RG_map_fw[t_qubit][l] == RED) {
                            // O qubit está num estado de superposição (2), então, não altere para 0 ou 1 sem verificar
                            RG_map_fw[t_qubit][l+1] = RED;
                        } else {
                            // Normalmente aplica a lógica da CNOT
                            RG_map_fw[t_qubit][l+1] = output[1];
                        }
                        
                        
                        //ver se o target e control ficam os 2 a vermelho
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_fw[nq][l+1]==RED){
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

/*printf("FORWARD COLOURING (RG_map_fw):\n");
	for (int nq = 0; nq < NQ; nq++) {
		printf("Qubit %d: ", nq);
		for (int l = 0; l <= L; l++) {
			printf("%d ", RG_map_fw[nq][l]);
		}
		printf("\n");
	}*/


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
                        colourT input_RG;
                        colourT output_RG;
                        
                        //get colors from RG_map_bw
                        //ver n�mero que est� no mapa para este qubit nesta layer: 0,1 ou 2
                        output_RG=RG_map_bw[qubit][l+1];
                        
                        //se for vermelho, sai
                        if (output_RG==RED) {
                            continue;
                        }
                        
                        //compute output colour for this gate
                        input_RG = gate_output_1(output_RG, name);
                        
                        // set colors in RG_map_bw
                        RG_map_bw[qubit][l]=input_RG;
                        
                        
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_bw[nq][l+1]==RED){
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
                        colourT input_RG;
                        colourT output_RG;
                        //get colors from RG_map_bw
                        //ver n�mero que est� no mapa para este qubit nesta layer: 0,1 ou 2
                        output_RG=RG_map_bw[qubit][l+1];
                        
                        //se for vermelho, sai
                        if (output_RG==RED) {
                            continue;
                        }
                        
                        //compute output colour for this gate
                        input_RG = gate_output_1(output_RG, name);
                        
                        // set colors in RG_map_bw
                        RG_map_bw[qubit][l]=input_RG;
                        
                        
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_bw[nq][l+1]==RED){
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
                        colourT output_RG_c = RG_map_bw[c_qubit][l + 1];
                        colourT output_RG_t = RG_map_bw[t_qubit][l + 1];
                        
                        // Se ambos forem 2 (irrelevantes), continua
                        if (output_RG_c == RED && output_RG_t == RED) {
                            continue;
                        }
                        
                        colourT output_RG[2] = {output_RG_c, output_RG_t};
                        colourT input[2];
                        gate_output_2(output_RG, name, input);  // Calcula input RG a partir de output RG
                        
                        // Atualiza RG_map_bw (estamos no backward sweep!)
                        RG_map_bw[c_qubit][l] = input[0];
                        
                        if (output_RG_t == RED) {
                            RG_map_bw[t_qubit][l] = RED;
                        } else {
                            RG_map_bw[t_qubit][l] = input[1];
                        }
                        
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_fw[nq][l+1]==RED){
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
                        colourT output_RG_c = RG_map_bw[c_qubit][l + 1];
                        colourT output_RG_t = RG_map_bw[t_qubit][l + 1];
                        
                        // Se ambos forem irrelevantes (2), não há nada a fazer
                        if (output_RG_c == RED && output_RG_t == RED) {
                            continue;
                        }
                        
                        colourT output_RG[2] = {output_RG_c, output_RG_t};
                        colourT input[2];
                        gate_output_2(output_RG, name, input);  // Calcula input RG a partir de output RG
                        
                        // Atualiza RG_map_bw (backward sweep!)
                        RG_map_bw[c_qubit][l] = input[0];
                        
                        if (output_RG_t == RED) {
                            RG_map_bw[t_qubit][l] = RED;
                        } else {
                            RG_map_bw[t_qubit][l] = input[1];
                        }
                        
                        for (int nq=0; nq<NQ;nq++){
                            if(RG_map_fw[nq][l+1]==RED){
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

/*printf("BACKWARD COLOURING (RG_map_bw):\n");
	for (int nq = 0; nq < NQ; nq++) {
		printf("Qubit %d: ", nq);
		for (int l = 0; l <= L; l++) {
			printf("%d ", RG_map_bw[nq][l]);
		}
		printf("\n");
	}*/
         
	//ver a cor e binary digit final
	
	for (int a=0; a<=L; a++){
		for (int b=0; b<NQ; b++){
            if (RG_map_fw[b][a]==RED && RG_map_bw[b][a]==RED){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
                RG_map[b][a]=RED;
			}
			
            else if (RG_map_fw[b][a]==RED){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=RG_map_bw[b][a];
			}
            else if (RG_map_bw[b][a]==RED){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=RG_map_fw[b][a];
			}
			
			else if (RG_map_fw[b][a]==RG_map_bw[b][a]){
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
				RG_map[b][a]=RG_map_fw[b][a];
			} else {
				//printf("fw[%d][%d]=%d, bw[%d][%d]=%d\n", b, a, RG_map_fw[b][a], b, a, RG_map_bw[b][a]);
                RG_map[b][a]=INVALID; //caso de conflito
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
    colourT* final_states_output = new colourT[(L-1)*NQ];
	//guardar no formato de um unico array

	for (int q = 0; q < NQ; q++) {
		for (int l = 1; l < L; l++) {
			final_states_output[(l - 1)*NQ+q] = RG_map[q][l];
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
