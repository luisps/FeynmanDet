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
	const int green=0;
	const int red=1;
 //nomear as branching gates de acordo com tipo GP
    int namesBG1P0[1] = {1}; //H
    int namesBG1P1[3] = {11, 12}; //rx,ry,rz
    int namesBG2P0[1] = {21}; //CX
    int *namesB[3] = {namesBG1P0,namesBG1P1,namesBG2P0};
    
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits; 
    const int N = (int)(powf(2.f, (float)NQ));
 //definir array com as cores para todos qubits e layers
    int coloursf[L][NQ]; //matriz com todos os lugares no estados interm�dios (forward sweep)
    int coloursb[L][NQ]; //matriz com todos os lugares no estados interm�dios (backward sweep)
    
    for (int q=0; q<NQ;q++){
    	for (int l=0; l<L; l++){
    		coloursf[l][q]=green; //mete tudo a 0, tudo a green
    		coloursb[l][q]=green;
		}
	}    
	//guardar a cor anterior
    int previousf[NQ];
	int previousb[NQ];	
		
	for (int i = 0; i < NQ; ++i) {
        previousf[i] = coloursf[0][i];
        previousb[i] = coloursb[L-1][i];
    }
    
    int inter_states[L-1][NQ];
    int inter_states_b[L-1][NQ];
    int previous_state;
    
    //FS
    for (int l=0; l<L;l++){ //FS: percorrer as layers - vai layer a layer
    	TCircuitLayer* layer = &circuit->layers[l];
    	
    	
        for (int GT = 0; GT < 4; GT++){ //percorrer os tipos de gates dentro de cada layer
        	switch (GT) {
				case 0: {
					TGate1P0 * g1p0_ptr = (TGate1P0 *) layer->gates[0];
					for (int g = 0; g < layer->num_type_gates[GT]; g++, g1p0_ptr++){
						int name = g1p0_ptr->name;
						int b=0;
						if (name==namesB[0][0]){
							b =1;
						}	
						//b fica a 1 se a gate for branching
						//FS: colorir de acordo o qubit anterior e natureza da gate
						//se for verde antes e a gate for n�o branching fica verde
						//se n�o fica vermelho
						if (previousf[g1p0_ptr->qubit]==green && b==0){
							coloursf[l+1][g1p0_ptr->qubit]=green;
						} else {coloursf[l+1][g1p0_ptr->qubit]=red;
						}
						//se for verde vamos ver o n�mero que tem 1 ou 0
						//e dizer o proximo numero conforme
						
						//avanca o estado
						if(l==0){
    						previous_state=init_state[g1p0_ptr->qubit];
						} else { 
							previous_state=inter_states[l-1][g1p0_ptr->qubit];
						}
						//se for vermelho, estado interm�dio fica vermelho
						if (coloursf[l+1][g1p0_ptr->qubit]==red){
							inter_states[l][g1p0_ptr->qubit]=2;
						//se n�o, vai-se ver a gate output	
						} else { inter_states[l][g1p0_ptr->qubit]=gate_output_1(previous_state,name);
						}	
					}
					break;
				}		
				//igual para os outros 3 casos
    			case 1:{
    				TGate1P1 * g1p1_ptr = (TGate1P1 *) layer->gates[1];
					for (int g = 0; g < layer->num_type_gates[GT]; g++, g1p1_ptr++){
						int name = g1p1_ptr->fdata.name;
						int b=0;
						for (int i=0;i<3;i++){
							if (name==namesB[1][i]){
								b=1;
							}	
						} 
						
						
						//FS: colorir de acordo o qubit anterior e natureza da gate
						if (previousf[g1p1_ptr->fdata.qubit]==green && b==0){
							coloursf[l+1][g1p1_ptr->fdata.qubit]=green;
						} else {coloursf[l+1][g1p1_ptr->fdata.qubit]=red;
						}
				
						if(l==0){
    						previous_state=init_state[g1p1_ptr->fdata.qubit];
    		
						} else { 
							previous_state=inter_states[l-1][g1p1_ptr->fdata.qubit];
						
						}
					
						
						if (coloursf[l+1][g1p1_ptr->fdata.qubit]==red){
							inter_states[l][g1p1_ptr->fdata.qubit]=2;
						
							
						} else { inter_states[l][g1p1_ptr->fdata.qubit]=gate_output_1(previous_state,name);
						
						}
				
			
					}    			
					break;
				}
	
     			case 2:{
     				TGate2P0 * g2p0_ptr = (TGate2P0 *) layer->gates[2];
     				for (int g = 0; g < layer->num_type_gates[GT]; g++, g2p0_ptr++){
     					int name = g2p0_ptr->name;
						int b=0;
						if (name==namesB[2][0]){
							b=1;
						}
						//FS: colorir de acordo o qubit anterior e natureza da gate
						if (previousf[g2p0_ptr->c_qubit]==green && previousf[g2p0_ptr->t_qubit]==green){
							coloursf[l+1][g2p0_ptr->c_qubit]=green;
							coloursf[l+1][g2p0_ptr->t_qubit]=green;
						} else if (previousf[g2p0_ptr->c_qubit]==green && previousf[g2p0_ptr->t_qubit]==red){
							coloursf[l+1][g2p0_ptr->c_qubit]=green;
							coloursf[l+1][g2p0_ptr->t_qubit]=red;
						} else {
							coloursf[l+1][g2p0_ptr->c_qubit]=red;
							coloursf[l+1][g2p0_ptr->t_qubit]=red;
						}
						
						
						int previous_state_c;
						int previous_state_t;
						
						if(l==0){
    						previous_state_c=init_state[g2p0_ptr->c_qubit];
    						previous_state_t=init_state[g2p0_ptr->t_qubit];
						} else { 
							previous_state_c=inter_states[l-1][g2p0_ptr->c_qubit];
							previous_state_t=inter_states[l-1][g2p0_ptr->t_qubit];
						}
						
						
						int previous_aux[2];
						previous_aux[0]=previous_state_c;
						previous_aux[1]=previous_state_t;
						
						
						if (coloursf[l+1][g2p0_ptr->c_qubit]==red){
							inter_states[l][g2p0_ptr->c_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states[l][g2p0_ptr->c_qubit]=output[0];
						}
					
						if (coloursf[l+1][g2p0_ptr->t_qubit]==red){
							inter_states[l][g2p0_ptr->t_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states[l][g2p0_ptr->t_qubit]=output[1];
						}
					
					
						
						//printf("This gates name is %d, and it is %d \n", name, b);
					}
					break;
				 }
     				
    			case 3:{
    				TGate2P1 * g2p1_ptr = (TGate2P1 *) layer->gates[3];
					for (int g = 0; g < layer->num_type_gates[GT]; g++, g2p1_ptr++){
						int b=0;
						int name=g2p1_ptr->fdata.name;
					
						//FS: colorir de acordo o qubit anterior e natureza da gate
						if (previousf[g2p1_ptr->fdata.c_qubit]==green && previousf[g2p1_ptr->fdata.t_qubit]==green){
							coloursf[l+1][g2p1_ptr->fdata.c_qubit]=green;
							coloursf[l+1][g2p1_ptr->fdata.t_qubit]=green;
						} else if (previousf[g2p1_ptr->fdata.c_qubit]==green && previousf[g2p1_ptr->fdata.t_qubit]==red){
							coloursf[l+1][g2p1_ptr->fdata.c_qubit]=green;
							coloursf[l+1][g2p1_ptr->fdata.t_qubit]=red;
						} else {
							coloursf[l+1][g2p1_ptr->fdata.c_qubit]=red;
							coloursf[l+1][g2p1_ptr->fdata.t_qubit]=red;
						}
						
						int previous_state_t;
						int previous_state_c;
						
						if(l==0){
    						previous_state_c=init_state[g2p1_ptr->fdata.c_qubit];
						} else { 
							previous_state_c=inter_states[l-1][g2p1_ptr->fdata.c_qubit];
						}
						
						
						
						if(l==0){
    						previous_state_t=init_state[g2p1_ptr->fdata.t_qubit];
						} else { 
							previous_state_t=inter_states[l-1][g2p1_ptr->fdata.t_qubit];
						}
						
						
						int previous_aux[2];
						previous_aux[0]=previous_state_c;
						previous_aux[1]=previous_state_t;
						
						if (coloursf[l+1][g2p1_ptr->fdata.c_qubit]==red){
							inter_states[l][g2p1_ptr->fdata.c_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states[l][g2p1_ptr->fdata.c_qubit]=output[0];
						}
					
						if (coloursf[l+1][g2p1_ptr->fdata.t_qubit]==red){
							inter_states[l][g2p1_ptr->fdata.t_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states[l][g2p1_ptr->fdata.t_qubit]=output[1];
						}
						
					
					}
					    				
					break;
				}
    			default:
        			printf("Gate type invalid");
			} 
        }
        // guardar o array de cores no previous p/ pode ser usado no loop seguinte
        for (int i = 0; i < NQ; i++) {
        	previousf[i] = coloursf[l+1][i];
        	
    	}
    }
    

	//BS
    for (int l=0; l<L;l++){ //percorrer as layers
    	TCircuitLayer* layer = &circuit->layers[L-(l+1)];
        for (int GT = 0; GT < 4; GT++){ //percorrer os tipos de gates
        	switch (GT) {
				case 0: {
					TGate1P0 * g1p0_ptr = (TGate1P0 *) layer->gates[0];
					for (int g = 0; g < layer->num_type_gates[GT]; g++, g1p0_ptr++){
						int name = g1p0_ptr->name;
						int b=0;
						if (name==namesB[0][0]){
							b =1;
						}	
						//BS:colorir de acordo o qubit anterior e natureza da gate
						if (previousb[g1p0_ptr->qubit]==green && b==0){
							coloursb[L-(l+2)][g1p0_ptr->qubit]=green;
						} else {coloursb[L-(l+2)][g1p0_ptr->qubit]=red;
						}
						
						
						//states
						if(l==0){
    						previous_state=final_state[g1p0_ptr->qubit];
						} else { 
							previous_state=inter_states_b[L-(l+1)][g1p0_ptr->qubit];
						}
						
						if (coloursb[L-(l+2)][g1p0_ptr->qubit]==red){
							inter_states_b[L-(l+2)][g1p0_ptr->qubit]=2;
							
						} else { inter_states_b[L-(l+2)][g1p0_ptr->qubit]=gate_output_1(previous_state,name);
						}
					}
				
					break;
				}		
    			case 1:{
    				TGate1P1 * g1p1_ptr = (TGate1P1 *) layer->gates[1];
					for (int g = 0; g < layer->num_type_gates[GT]; g++, g1p1_ptr++){
						int name = g1p1_ptr->fdata.name;
						int b=0;
						for (int i=0;i<3;i++){
							if (name==namesB[1][i]){
								b=1;
							}	
						} 
					
						//BS:colorir de acordo o qubit anterior e natureza da gate
						if (previousb[g1p1_ptr->fdata.qubit]==green && b==0){
							coloursb[L-(l+2)][g1p1_ptr->fdata.qubit]=green;
						} else {coloursb[L-(l+2)][g1p1_ptr->fdata.qubit]=red;
						}
						
						
						if(l==0){
    						previous_state=final_state[g1p1_ptr->fdata.qubit];
    		
						} else { 
							previous_state=inter_states_b[L-(l+1)][g1p1_ptr->fdata.qubit];
						
						}
					
						
						if (coloursb[L-(l+2)][g1p1_ptr->fdata.qubit]==red){
							inter_states_b[L-(l+2)][g1p1_ptr->fdata.qubit]=2;
						
							
						} else { inter_states_b[L-(l+2)][g1p1_ptr->fdata.qubit]=gate_output_1(previous_state,name);
						
						}
						
					}    			
					break;
				}
	
     			case 2:{
     				TGate2P0 * g2p0_ptr = (TGate2P0 *) layer->gates[2];
     				for (int g = 0; g < layer->num_type_gates[GT]; g++, g2p0_ptr++){
     					int name = g2p0_ptr->name;
						int b=0;
						if (name==namesB[2][0]){
							b=1;
						}
					
						//BS:colorir de acordo o qubit anterior e natureza da gate
						if (previousb[g2p0_ptr->c_qubit]==green && previousb[g2p0_ptr->t_qubit]==green){
							coloursb[L-(l+2)][g2p0_ptr->c_qubit]=green;
							coloursb[L-(l+2)][g2p0_ptr->t_qubit]=green;
						} else if (previousb[g2p0_ptr->c_qubit]==green && previousb[g2p0_ptr->t_qubit]==red){
							coloursb[L-(l+2)][g2p0_ptr->c_qubit]=green;
							coloursb[L-(l+2)][g2p0_ptr->t_qubit]=red;
						} else {
							coloursb[L-(l+2)][g2p0_ptr->c_qubit]=red;
							coloursb[L-(l+2)][g2p0_ptr->t_qubit]=red;
						}
						
						
						
						int previous_state_c;
						int previous_state_t;
						
						if(l==0){
    						previous_state_c=final_state[g2p0_ptr->c_qubit];
    						previous_state_t=final_state[g2p0_ptr->t_qubit];
						} else { 
							previous_state_c=inter_states_b[L-(l+1)][g2p0_ptr->c_qubit];
							previous_state_t=inter_states_b[L-(l+1)][g2p0_ptr->t_qubit];
						}
						
						
						int previous_aux[2];
						previous_aux[0]=previous_state_c;
						previous_aux[1]=previous_state_t;
						
						
						if (coloursb[L-(l+2)][g2p0_ptr->c_qubit]==red){
							inter_states_b[L-(l+2)][g2p0_ptr->c_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states_b[L-(l+2)][g2p0_ptr->c_qubit]=output[0];
						}
					
						if (coloursb[L-(l+2)][g2p0_ptr->t_qubit]==red){
							inter_states_b[L-(l+2)][g2p0_ptr->t_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states_b[L-(l+2)][g2p0_ptr->t_qubit]=output[1];
						}
					
						
						//printf("This gates name is %d, and it is %d \n", name, b);
					}
					break;
				 }
     				
    			case 3:{
    				TGate2P1 * g2p1_ptr = (TGate2P1 *) layer->gates[3];
					for (int g = 0; g < layer->num_type_gates[GT]; g++, g2p1_ptr++){
						int b=0;
						int name=g2p1_ptr->fdata.name;
						
						//BS:colorir de acordo o qubit anterior e natureza da gate
						if (previousb[g2p1_ptr->fdata.c_qubit]==green && previousb[g2p1_ptr->fdata.t_qubit]==green){
							coloursb[L-(l+2)][g2p1_ptr->fdata.c_qubit]=green;
							coloursb[L-(l+2)][g2p1_ptr->fdata.t_qubit]=green;
						} else if (previousb[g2p1_ptr->fdata.c_qubit]==green && previousb[g2p1_ptr->fdata.t_qubit]==red){
							coloursb[L-(l+2)][g2p1_ptr->fdata.c_qubit]=green;
							coloursb[L-(l+2)][g2p1_ptr->fdata.t_qubit]=red;
						} else {
							coloursb[L-(l+2)][g2p1_ptr->fdata.c_qubit]=red;
							coloursb[L-(l+2)][g2p1_ptr->fdata.t_qubit]=red;
						}
						
						
						int previous_state_t;
						int previous_state_c;
						
						if(l==0){
    						previous_state_c=final_state[g2p1_ptr->fdata.c_qubit];
    						previous_state_t=final_state[g2p1_ptr->fdata.t_qubit];
						} else { 
							previous_state_c=inter_states_b[L-(l+1)][g2p1_ptr->fdata.c_qubit];
							previous_state_t=inter_states_b[L-(l+1)][g2p1_ptr->fdata.t_qubit];
						}
						
						
						int previous_aux[2];
						previous_aux[0]=previous_state_c;
						previous_aux[1]=previous_state_t;
						
						if (coloursb[L-(l+2)][g2p1_ptr->fdata.c_qubit]==red){
							inter_states_b[L-(l+2)][g2p1_ptr->fdata.c_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states_b[L-(l+2)][g2p1_ptr->fdata.c_qubit]=output[0];
						}
					
						if (coloursb[L-(l+2)][g2p1_ptr->fdata.t_qubit]==red){
							inter_states_b[L-(l+2)][g2p1_ptr->fdata.t_qubit]=2;
							
						} else { 
							int output[2];
							gate_output_2(previous_aux,name, output);
							inter_states_b[L-(l+2)][g2p1_ptr->fdata.t_qubit]=output[1];
						}
						
						
					}
					    				
					break;
				}
    			default:
        			printf("Gate type invalid");
			} 
        }
        //defining previous colours 
        for (int i = 0; i < NQ; i++) {
        	previousb[i] = coloursb[L-(l+2)][i];
    	}
    }
    
    /* printf("This is forward sweep \n");
    for (int q=0; q<NQ;q++){
    	for (int l=0; l<L; l++){
    		printf("%d ", coloursf[l][q]);
		}
		printf("\n");
	}
	
	//printf("This is backwards sweep \n");
    for (int q=0; q<NQ;q++){
    	for (int l=0; l<L; l++){
    		printf("%d ", coloursb[l][q]);
		}
		printf("\n");
	}
	*/
	int final_colour[L-1][NQ];
	
	int* final_colour_output = new int[(L-1)*NQ];
	
	//printf("This is the full sweep \n");
	
	//ver a cor final (um vez o outro, pois 1x0 d� 0: ou seja, se for green num fica green)
	for (int b=0; b<NQ; b++){
		for (int a=0; a<L-1; a++){
			final_colour_output[a*(L-1)+b] = coloursf[a+1][b] * coloursb[a][b];	
			final_colour[a][b] = coloursf[a+1][b] * coloursb[a][b];
			//printf("%d ", final_colour[a][b]);
			
		}
		//printf("\n");
	
	}
	
	/*
	printf("This is intermediate state forward \n");
    for (int q=0; q<NQ;q++){
    	for (int l=0; l<L-1; l++){
    		printf("%d ", inter_states[l][q]);
		}
		printf("\n");
	}
	
	
	printf("This is intermediate state backward \n");
	for (int q=0; q<NQ;q++){
    	for (int l=0; l<L-1; l++){
    		printf("%d ", inter_states_b[l][q]);
		}
		printf("\n");
	}
	*/
	//atribuir aos estados verdes o seu valor (1 ou 0)
	// se for verdes mais diferentes � imposs�vel
	int final_states[L-1][NQ];
	
	for(int l=0; l<L-1;l++){
		for(int q=0; q<NQ; q++){
			
			if(coloursf[l+1][q]==green){
				final_states[l][q]=inter_states[l][q];
			}
			if(coloursb[l][q]==green){
				final_states[l][q]=inter_states_b[l][q];
			}
			if(coloursb[l][q]==green && coloursf[l+1][q]==green && inter_states[l][q]!=inter_states_b[l][q]){
				final_states[l][q]=-1;
			}
			if(final_colour[l][q]==red){
				final_states[l][q]=2;
			}
		

		}
	}
	
	int* final_states_output = new int[(L-1)*NQ];
	//guardar no formato de um �nico array
	
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








