//
//  simulate_useful.cpp
//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
//FUNCIONA MAS � ANTIGA (SEM START LAYER, ETC)


#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include "layer.hpp"
#include "layer.cpp"
#include "simulation_attempt.cpp"


void simulate_useful_paths (TCircuit *circuit,
                         int init_state, int final_state,
                         float& aR, float& aI) {
    
    
    const int L = circuit->size->num_layers;
    int *ndxs = new int[L-1];
    const int N = (int)(powf(2.f, (float)circuit->size->num_qubits));
    const int NQ=circuit->size->num_qubits;


	int* states = new int[(L-1)*NQ];
	int init_state_arr[NQ];
	int final_state_arr[NQ];
	
	for (int i=0; i<NQ; i++){
		init_state_arr[i]=qb_value(i,init_state);
		final_state_arr[i]=qb_value(i,final_state);	
    }
	
	states=fs_bs_RG(circuit, init_state_arr, final_state_arr);
    
    float sumR=0.f, sumI=0.f;
    
    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;
    
    //printf("These are the states: \n");

	for(int s=0; s<(L-1)*NQ; s++){
        //printf("%d ", states[s]);
    	if(states[s]==-1){
    			states[0]=-1;
			break;
		}
	}

	//printf("\n");
    if(states[0]==-1){
    	printf("IMPOSSIBLE: EXITING IMMEDIATELY\n");
    	aR=0.f;
    	aI=0.f;
    	printf ("< %d | U | %d > = %.6f + i %.6f\n", final_state, init_state, aR, aI);
    	return;
	} else{
	

    while (ndxs[0] < N) {
    	
		// printf("%d ", ndxs[0]);        
        float pathR=1.f;
		float pathI=0.f;
                
        int current_state = init_state;
        int next_state;
        //int next_state_cl[NQ];
        //int next_state_arr[NQ];
        
        for (int l=0 ; l<L ; l++) {
            float lR=1.f;
			float lI=0.f;
			int is_zero = 0;
            
            next_state = (l< L-1 ? ndxs[l] : final_state);
            
            //printf("%d ", next_state);
            if(l < L-1){
            	for (int i=0; i<NQ; i++){
            		const int g_vs_r = states[i*(L-1)+l];
					const int  qb =qb_value(i,next_state); //guarda ndxs em binario
					//vai meter aqui o valor do fs_bs_RG
					//next_state_cl[i]=states[i*(L-1)+l]; //guarda o que tem de ser de acordo com a funcao
					//se nao for vermelho e forem diferentes, caminho da 0
					//se tivermos um verde que tem de ser 1 mas ha um ndxs que tem la 0
						//deita-se esse logo fora

                    //next_state_arr[i]=qb_value(i,next_state); 
					//next_state_cl[i]=states[i*(L-1)+l]; 
					//if(next_state_cl[i]!=2 && next_state_cl[i]!=next_state_arr[i]){

					if(g_vs_r!=2 && g_vs_r!=qb){
						pathR = pathI = 0.f;
                        //printf("States were different! -> layer %d\n", l);
						is_zero=1;
						break;
					}
    			}
    		}
    		
    		if(is_zero==1){
    			break;;
			}
    		
            //printf("States match! -> layer %d\n", l);

            TCircuitLayer *layer = &circuit->layers[l];

            layer_w(layer, l, current_state, next_state, lR, lI);
        	
            if (complex_abs_square(lR, lI) <= 0.f) {
                pathR = pathI = 0.f;
                break;
            }
            complex_multiply(pathR, pathI, lR, lI, pathR, pathI);
            current_state = next_state;

        }
        //printf("\n ");
        //printf("%d + i %d \n", pathR, pathI);
        sumR += pathR;
        sumI += pathI;

        
        int ll;
        for (ll=L-2 ; ll>=0 ; ll--) {
            ndxs[ll]++;
            if (ndxs[ll]==N && ll>0)  {
                ndxs[ll] = 0;
            }
            else
                break;
        }

    }

    aR = sumR;
    aI = sumI;
    //printf ("< %d | U | %d > = %.6f + i %.6f\n", final_state, init_state, aR, aI);
    printf ("< %d | U | %d > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    }

}