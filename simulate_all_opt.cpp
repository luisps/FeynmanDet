//
//  simulate_useful.cpp
//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//

#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include "layer.hpp"
#include "layer.cpp"
#include "funcao.cpp"

#define NEW_VERSION

//#define PRINT_STATE

void simulate_all_paths (TCircuit *circuit,
                         int init_state, int final_state,
                         float& aR, float& aI) {
    
    //printf ("simulate <%d | U | %d>\n", init_state, final_state);
    //return;
    const int L = circuit->size->num_layers;
    int *ndxs = new int[L-1];
    const int N = (int)(powf(2.f, (float)circuit->size->num_qubits));
    const int NQ=circuit->size->num_qubits;
#ifdef NEW_VERSION       
    float* wR=new float[L-1]; //NEW
    float* wI=new float[L-1]; //NEW
#endif

	int* states = new int[(L-1)*NQ];
	int init_state_arr[NQ];
	int final_state_arr[NQ];
	
	for (int i=0; i<NQ; i++){
		init_state_arr[i]=qb_value(i,init_state);
		final_state_arr[i]=qb_value(i,final_state);	
    }
	

	//os estados s�o dados pela fun��o (com n�mero 1 0 2 ou -1)
    //states=fs_bs_RG(circuit, init_state_arr, final_state_arr);
   
    for (int q=0; q<NQ;q++){
    	for (int l=0; l<L-1; l++){
    		states[q*(L-1)+l]=2;
		}
	}

    float sumR=0.f, sumI=0.f;
#ifdef NEW_VERSION    
    int start_layer=0; //NEW
#endif
    //inicializar estados interm�dios
    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;
    
    //se for -1 em algum lado mete -1 no in�cio
	for(int s=0; s<(L-1)*NQ; s++){
    	//printf("%d ", states[s]);
    	if(states[s]==-1){
    		printf("IMPOSSIBLE: EXITING IMMEDIATELY\n");
    		aR=0.f;
    		aI=0.f;
    		printf ("< %d | U | %d > = %.6f + i %.6f\n", final_state, init_state, aR, aI);
    		return;
		}
	}

	//para todos os estados interm�dios
    while (ndxs[0] < N) {
    	
		// printf("%d ", ndxs[0]);
        // .... process path : init_state-> ndxs[0] -> ..... -> ndxs[L-2] -> final_state
        int next_state,l;
        
#ifdef NEW_VERSION        
       /* for(l=0; l<start_layer; l++) {
        	
			complex_multiply (pathR, pathI, wR[l], wI[l], pathR, pathI);
		}*/
		float pathR = (start_layer==0? 1.f : wR[start_layer-1]);
		float pathI = (start_layer==0? 0.f : wI[start_layer-1]);
        int current_state = (start_layer==0? init_state : ndxs[start_layer-1]);
#else
        float pathR=1.f, pathI=0.f;
        int current_state = init_state;
#endif        
        //int next_state_cl[NQ];
        //int next_state_arr[NQ];
        
        bool zero_weight_layer=false; 
        
        
#ifdef NEW_VERSION        
        for (l=start_layer ; l<L ; l++) {
#else
        for (l=0 ; l<L ; l++) {
#endif        
            float lR=1.f, lI=0.f;
			int is_zero = 0;
            
            next_state = (l< L-1 ? ndxs[l] : final_state);
            
            //printf("%d ", next_state);
            if(l < L-1){
            	for (int i=0; i<NQ; i++){
            		const int g_vs_r = states[l*NQ +i];
					const int  qb =qb_value(i,next_state); //guarda ndxs em bin�rio
					//vai meter aqui o valor do fs_bs_RG
					//next_state_cl[i]=states[i*(L-1)+l]; //guarda o que tem de ser de acordo com a fun��o
					//se n�o for vermelho e forem diferentes, caminho d� 0
					//se tivermos um verde que tem de ser 1 mas h� um ndxs que tem l� 0
						//deita-se esse logo fora
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
    			break;
			}
    		
			//printf("States match! -> layer %d\n", l);
    		
    		//calcular probabilidades
            TCircuitLayer *layer = &circuit->layers[l];
            //layer_w atualiza valor de lR e lI
            layer_w(layer, l, current_state, next_state, lR, lI);
        	
        	//faz mult e avan�a para o proximo estado
            complex_multiply(pathR, pathI, lR, lI, pathR, pathI);
            
#ifdef NEW_VERSION       
            wR[l]=pathR;
            wI[l]=pathI;
#endif     
        	//se prob for negativa, est� mal
            if (complex_abs_square(lR, lI) <= 0.f) {
            	zero_weight_layer=true;
                pathR = pathI = 0.f;
                break;
            }
         
			current_state = next_state;	
        }
        //printf("\n ");
        //printf("%d + i %d \n", pathR, pathI);
        sumR += pathR;
        sumI += pathI;

        
        // iterate to next path: mudar o ndxs
        int ll;
        for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=0 ; ll--) {
            ndxs[ll]++;
#ifdef NEW_VERSION        
            start_layer=ll;
#endif            
            if (ndxs[ll]==N && ll>0)  {  // this loop resets
                ndxs[ll] = 0;
            }
            else
                break;
        }

    }

    aR = sumR;
    aI = sumI;
    printf ("< %d | U | %d > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
#ifdef PRINT_STATE        
    printf ("< %d | U | %d > = %.6f + i %.6f\n", final_state, init_state, aR, aI);
#endif
    
}
