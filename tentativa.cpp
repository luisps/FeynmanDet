#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include "layer.hpp"
#include "circuit.cpp"

//int forwardSweepRG(TCircuit *circuit){
 //   const int L = circuit->size->num_layers;
   // const int NQ = circuit->size->num_qubits; 
  //  const int N = (int)(powf(2.f, (float)NQ));   
 //   int colours[L-1][NQ];
    //percorrer layer a layer
   // for (int l=0; l<L-1;l++){
    //    TCircuitLayer *layer = &circuit->layers[l];
    //    int ng = layer -> num_gates;
      //  for (int g=0; g<ng; g++){
            //if(gate[0]=='S' or gate[0]=='T' or gate[0]=='I')
      //      if()
      //  }
  //  }


int forwardSweepRG(TCircuit *circuit){
	const int green=0;
	const int red=1;
 //nomear as branching gates de acordo com tipo GP
    int namesBG1P0[1] = {1};
    int namesBG1P1[3] = {11, 12, 13};
    int namesBG2P0[1] = {21};
    int *namesB[3] = {namesBG1P0,namesBG1P1,namesBG2P0};
    
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits; 
    const int N = (int)(powf(2.f, (float)NQ));
 //definir array com as cores para todos qubits e layers
    int colours[L][NQ]; //matriz com todos os lugares no estados intermédios
    //debugging
    printf("This is colours before everything \n");
    for (int q=0; q<NQ;q++){
    	for (int l=0; l<L; l++){
    		colours[l][q]=0;
    		printf("%d ", colours[l][q]);
		}
		printf("\n");
	}    
	
	
    int previous[NQ];		
	for (int i = 0; i < NQ; ++i) {
        previous[i] = colours[0][i];
    }
    
    for (int l=0; l<L;l++){ //percorrer as layers
    	TCircuitLayer* layer = &circuit->layers[l];
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
						if (previous[g1p0_ptr->qubit]==green && b==0){
							colours[l+1][g1p0_ptr->qubit]=green;
						} else {colours[l+1][g1p0_ptr->qubit]=red;
						}
						printf("This gates name is %d, and it is %d \n", name, b);
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
						if (previous[g1p1_ptr->fdata.qubit]==green && b==0){
							colours[l+1][g1p1_ptr->fdata.qubit]=green;
						} else {colours[l+1][g1p1_ptr->fdata.qubit]=red;
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
						if (previous[g2p0_ptr->c_qubit]==green && previous[g2p0_ptr->t_qubit]==green){
							colours[l+1][g2p0_ptr->c_qubit]=green;
							colours[l+1][g2p0_ptr->t_qubit]=green;
						} else if (previous[g2p0_ptr->c_qubit]==green && previous[g2p0_ptr->t_qubit]==red){
							colours[l+1][g2p0_ptr->c_qubit]=green;
							colours[l+1][g2p0_ptr->t_qubit]=red;
						} else {
							colours[l+1][g2p0_ptr->c_qubit]=red;
							colours[l+1][g2p0_ptr->t_qubit]=red;
						}
						
						printf("This gates name is %d, and it is %d \n", name, b);
					}
					break;
				 }
     				
    			case 3:{
    				TGate2P1 * g2p1_ptr = (TGate2P1 *) layer->gates[3];
					for (int g = 0; g < layer->num_type_gates[GT]; g++, g2p1_ptr++){
						int b=0;
						if (previous[g2p1_ptr->fdata.c_qubit]==green && previous[g2p1_ptr->fdata.t_qubit]==green){
							colours[l+1][g2p1_ptr->fdata.c_qubit]=0;
							colours[l+1][g2p1_ptr->fdata.t_qubit]=0;
						} else if (previous[g2p1_ptr->fdata.c_qubit]==green && previous[g2p1_ptr->fdata.t_qubit]==red){
							colours[l+1][g2p1_ptr->fdata.c_qubit]=0;
							colours[l+1][g2p1_ptr->fdata.t_qubit]=1;
						} else {
							colours[l+1][g2p1_ptr->fdata.c_qubit]=1;
							colours[l+1][g2p1_ptr->fdata.t_qubit]=1;
						}
					}
					    				
					break;
				}
    			default:
        			printf("Gate type invalid");
			} 
        }
        for (int i = 0; i < NQ; i++) {
        	previous[i] = colours[l+1][i];
    	}
    }
    for (int q=0; q<NQ;q++){
    	for (int l=0; l<L; l++){
    		printf("%d ", colours[l][q]);
		}
		printf("\n");
	}
    return 0;
}


//single qubit gates:
    //estado anterior green -> NB gate -> estado a green
    //estado anterior red -> NB gate -> estado a red
    //estado anterior green/red -> B gate -> estado a red
// 2 qubit gates:
    //t+c green ->all green
    //t+c red -> all red
    //t red, c green -> t red, c green
    // t green, c red -> all red


    int main() {
    const char * fileName = "circuit_211.data"; 

    // Call the read_circuit function
    TCircuit * circuit = read_circuit(fileName);

    // Check if the function succeeded
    if (circuit != NULL) {
        print_circuit(circuit);
        printf("This is a test for the FS: \n");
        int colours = forwardSweepRG(circuit);

        // Free the allocated memory
        free(circuit->size);
        free(circuit->layers);
        free(circuit);
    } else {
        fprintf(stderr, "Error reading circuit from file.\n");
    }

    printf("This is a test for the binary conversion: \n");

    for (int i=10; i>=0; i--){
        printf("%d", qb_value(i,13));
    }
    printf("\n");
    system("pause");
    return 0;
}
