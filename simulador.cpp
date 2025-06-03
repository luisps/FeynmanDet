//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#include "simulate_all.hpp"
#include <math.h>
#include <algorithm>
#include "complex.h"
#include "layer.hpp"
#include "layer.cpp"
//#include "funcao.cpp"
#include "colouring_LPSsuggestion.cpp"
#include <time.h>
#include <iostream> 

#define COLOURING
#define STL

void simulate_useful_paths (TCircuit *circuit, int init_state, int final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    int *ndxs = new int[L-1];
    const int N = (int)(powf(2.f, (float)circuit->size->num_qubits));
    const int NQ=circuit->size->num_qubits;

    float* wR=new float[L-1];
    float* wI=new float[L-1];

    int* states = new int[(L-1)*NQ];

    int init_state_arr[NQ];
    int final_state_arr[NQ];
    float sumR=0.f, sumI=0.f;
    for (int i=0; i<NQ; i++){
        init_state_arr[i]=qb_value(i,init_state);
        final_state_arr[i]=qb_value(i,final_state);	
    }

#ifdef COLOURING    
    states=fs_bs_RG(circuit, init_state_arr, final_state_arr);
#else
    for (int q=0; q<NQ;q++){
        for (int l=0; l<L-1; l++){
            states[q*(L-1)+l]=2;
        }
    }
#endif

#ifdef STL
    int start_layer=0;
#endif

    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;
    for(int s=0; s<(L-1)*NQ; s++){
        if(states[s]==-1){
            aR=0.f;
            aI=0.f;
            printf ("< %d | U | %d > = %.6f + i %.6f\n", final_state, init_state, aR, aI);
            return;
        }
    }

    while (ndxs[0] < N) {

        bool skip_path = false;

#ifdef STL 
        float pathR = (start_layer==0? 1.f : wR[start_layer-1]);
        float pathI = (start_layer==0? 0.f : wI[start_layer-1]);
        int current_state = (start_layer==0? init_state : ndxs[start_layer-1]);
#else
        float pathR=1.f; 
        float pathI=0.f;
        int current_state = init_state; 
#endif

#ifdef COLOURING
#ifdef STL
        if (start_layer > 0) {
            bool colouring_mismatch = false;
            for (int l_check = 0; l_check < start_layer; l_check++) {
                int expected_next_state = (l_check < L-1 ? ndxs[l_check] : final_state);
                for (int q = 0; q < NQ; q++) {
                    int expected = states[q*(L-1) + l_check];
                    if (expected == 2) continue; // No constraint on this qubit at this layer
                    int actual = qb_value(q, expected_next_state);
                    if (expected != actual) {
                        colouring_mismatch = true;
                        break;
                    }
                }
                if (colouring_mismatch) break;
            }
            if (colouring_mismatch) {
                // Path violates colouring before start_layer: skip immediately
                pathR = pathI = 0.f;
                skip_path = true;
            }
        }
#endif
#endif

        int next_state,l;
        int next_state_cl[NQ];
        int next_state_arr[NQ];
        bool zero_weight_layer=false;

#ifdef STL
        for (l=start_layer ; l<L ; l++) {
#else
        //for (l=0 ; l<L ; l++) { 
#endif
            float lR=1.f;
            float lI=0.f;
            int is_zero = 0;
            next_state = (l< L-1 ? ndxs[l] : final_state);

            if(l < L-1){
                for (int i=0; i<NQ; i++){
                    next_state_arr[i]=qb_value(i,next_state);
                    next_state_cl[i]=states[i*(L-1)+l];
                    if(next_state_cl[i]!=2 && next_state_cl[i]!=next_state_arr[i]){
                        pathR = pathI = 0.f;
                        is_zero=1;
                        break;
                    }
                }
            }
            if(is_zero==1){
                break;
            }
        
            TCircuitLayer *layer = &circuit->layers[l];
            layer_w(layer, l, current_state, next_state, lR, lI);
            complex_multiply(pathR, pathI, lR, lI, pathR, pathI);


#ifdef STL
            wR[l]=pathR;
            wI[l]=pathI;
#else
            wR[l]=lR;
            wI[l]=lI;
#endif
            if (complex_abs_square(lR, lI) <= 0.f) {
                zero_weight_layer=true;
                pathR = pathI = 0.f;
                break;
            }
            current_state = next_state;

        }
        if (!skip_path) {
            sumR += pathR;
            sumI += pathI;
        }

        int ll;
        for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=0 ; ll--) {
            ndxs[ll]++;
#ifdef STL
            start_layer=ll;
#endif
            if (ndxs[ll]==N && ll>0)  { 
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