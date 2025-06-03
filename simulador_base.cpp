//
//  simulate_useful.cpp
//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
                                                        //Colouring
#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include "layer.hpp"
#include "layer.cpp"
#include "funcao.cpp"
//#include "simulation_attempt_f.cpp"
//#include "simulation_attempt_03_10.cpp"
//#include "simulation_attempt_14_12.cpp"
//#include "colouring_LPSsuggestion.cpp"
//#include "simulation_attempt_funcional.cpp"
#include <time.h>




void simulate_useful_paths (TCircuit *circuit,
                         int init_state, int final_state,
                         float& aR, float& aI) {
    
    //queremos pegar num estado inicial e num estado final
    //calcular a parte cl�ssica com a fun��o
    //checkar se tem -1's -> descartar
    //s� acrescentar � soma se a parte cl�ssica dos estados interm�dios concordar com o que temos
    
    //clock_t start_pre, end_pre;
    
	//start_pre=clock();
	
    const int L = circuit->size->num_layers;
    int *ndxs = new int[L-1];
    const int N = (int)(powf(2.f, (float)circuit->size->num_qubits));
    const int NQ=circuit->size->num_qubits;

    float* wR=new float[L-1];//NEW
    float* wI=new float[L-1];//NEW

	int* states = new int[(L-1)*NQ];
	int init_state_arr[NQ];
	int final_state_arr[NQ];
	
    
    for (int i=0 ; i<L-1 ; i++) wR[i]=wI[i]=0.f; //indica que constante e float //NEW
	
	//qb_value guarda num array os estados final e inicial (que s�o n�meros decimais inicialmente)
	//qb_value de i,init_state=11 d� 1011
	for (int i=0; i<NQ; i++){
		init_state_arr[i]=qb_value(i,init_state);
		final_state_arr[i]=qb_value(i,final_state);	
    }
    
    
	//os estados sao dados pela funcao (com numero 1 0 2 ou -1)
	states=fs_bs_RG(circuit, init_state_arr, final_state_arr); //DIFF

    //end_pre=clock();
    //double time_pre=double(end_pre - start_pre)/double(CLOCKS_PER_SEC);
	//printf("Time taken pre-processing is: %lf \n", time_pre);
    
    
    //clock_t start_pos, end_pos;
    
	//start_pos=clock();
    
    float sumR=0.f, sumI=0.f;
    
    //inicializar estados intermedios
    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;
    
    
    //printf("These are the states: \n");
    
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
        
        float pathR=1.f; //amp de medir estados inicial tendo estado inicial=1
		float pathI=0.f;
                
        int current_state = init_state;
        int next_state,l;
        int next_state_cl[NQ];
        int next_state_arr[NQ];
        
        bool zero_weight_layer=false;
        
        
        for (l=0 ; l<L ; l++) { //percorre layers
            float lR=1.f;
			float lI=0.f;
			int is_zero = 0;
            
            next_state = (l< L-1 ? ndxs[l] : final_state);
            
            //printf("%d ", next_state);
            // ver quais s�o os estados possiveis e impossiveis e s� correr os que interessam
             if(l < L-1){
            	for (int i=0; i<NQ; i++){
					next_state_arr[i]=qb_value(i,next_state); //guarda ndxs em bin�rio
					//vai meter aqui o valor do fs_bs_RG
					next_state_cl[i]=states[i*(L-1)+l]; //guarda o que tem de ser de acordo com a fun��o
					//se n�o for vermelho e forem diferentes, caminho d� 0
					//se tivermos um verde que tem de ser 1 mas h� um ndxs que tem l� 0
						//deita-se esse logo fora
					if(next_state_cl[i]!=2 && next_state_cl[i]!=next_state_arr[i]){
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
            
        	wR[l]=lR;
            wI[l]=lI;
            
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

    }
    //end_pos=clock();
    
    //double time_pos=double(end_pos - start_pos)/double(CLOCKS_PER_SEC);
	// printf("Time taken running the circuit is: %lf \n", time_pos);