//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#include <math.h>
#include <algorithm>
#include "my_complex.h"
#include "layer.hpp"
#include "colouring_RG.hpp"
#include <time.h>
#include <iostream>

#include "simulator_RG.hpp"

static void print_RG (int *states, int NQ, int L) {
    fprintf (stderr, "Colouring: \n");
    for (int q = 0; q < NQ; q++) {
        fprintf (stderr, "qb%d: ", q);
        for (int l = 1; l < L; l++) {
            int colour= states[q * (L - 1) + (l - 1)];
            char s[5];
            switch (colour) {
                case -1:
                    snprintf (s,4,"I");
                    break;
                case 0:
                    snprintf (s,4,"G0");
                    break;
                case 1:
                    snprintf (s,4,"G1");
                    break;
                case 2:
                    snprintf (s,4,"R");
                    break;
            }
            fprintf (stderr, "%s\t", s);
        }
        fprintf (stderr, "\n");
    }
    fprintf (stderr, "\n");
}

void simulate_RG_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    StateT *ndxs = new StateT[L-1];
    const int NQ=circuit->size->num_qubits;
    const StateT N = 1 << NQ;
    StateT path_counter=0, path_NZ_counter=0;

    float* wR=new float[L-1];
    float* wI=new float[L-1];


    int init_state_arr[NQ];
    int final_state_arr[NQ];
    float sumR=0.f, sumI=0.f;
    for (int i=0; i<NQ; i++){
        init_state_arr[i]=qb_value(i,init_state);
        final_state_arr[i]=qb_value(i,final_state);	
    }

    // Do the colouring
    // 0 -> Green = 0
    // 1 -> Green = 1
    // 2 -> Red (branching)
    // -1 -> forward green != backward green : impossible

    // to store the colouring results
    int* colours;
    colours=fs_bs_RG(circuit, init_state_arr, final_state_arr);

    print_RG (colours, NQ, L);
    int start_layer=0;

    // verify if there is a -1 in the colouring
    // if yes, amplitude = 0 + 0j
    for(int s=0; s<(L-1)*NQ; s++){
        if(colours[s]==-1){
            aR=0.f;
            aI=0.f;
            printf ("EARLY TERMINATION: < %llu | U | %llu > = %.6f + i %.6f\n", (unsigned long long)final_state,
                (unsigned long long)init_state, aR, aI);
            return;
        }
    }

    // Simulation starts

    // all intermediate layers indexes to 0
    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;

    // main simulation loop
    while (ndxs[0] < N) {

        /*
        fprintf (stderr, "ndxs= ");
        for (int lll=0; lll<L-1 ; lll++)
            fprintf (stderr, "%llu ", ndxs[lll]);
        fprintf (stderr, "\n");*/
        
        bool skip_path = false;

        float pathR = (start_layer==0? 1.f : wR[start_layer-1]);
        float pathI = (start_layer==0? 0.f : wI[start_layer-1]);
        StateT current_state = (start_layer==0? init_state : ndxs[start_layer-1]);

        if (start_layer > 0) {
            bool colouring_mismatch = false;
            for (int l_check = 0; l_check < start_layer; l_check++) {
                StateT expected_next_state = (l_check < L-1 ? ndxs[l_check] : final_state);
                for (int q = 0; q < NQ; q++) {
                    int expected = colours[q*(L-1) + l_check];
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

        int l;
        StateT next_state;
        int next_state_cl[NQ];
        int next_state_arr[NQ];
        bool zero_weight_layer=false;
        bool is_zero=false;
        
        // iterate over layers
        for (l=start_layer ; l<L ; l++) {
            float lR=1.f;
            float lI=0.f;
            is_zero = false;
            next_state = (l< L-1 ? ndxs[l] : final_state);

            // verify whether this next state is allowed by colouring
            if(l < L-1){
                for (int i=0; i<NQ; i++){
                    next_state_arr[i]=qb_value(i,next_state);
                    next_state_cl[i]=colours[i*(L-1)+l];
                    if(next_state_cl[i]!=2 && next_state_cl[i]!=next_state_arr[i]){
                        pathR = pathI = 0.f;
                        is_zero=true;
                        break;
                    }
                }
            }
            if(is_zero){
                break;
            }
        
            TCircuitLayer *layer = &circuit->layers[l];
            layer_w(layer, l, current_state, next_state, lR, lI);
            complex_multiply(pathR, pathI, lR, lI, pathR, pathI);


            wR[l]=pathR;
            wI[l]=pathI;
            if (complex_abs_square(lR, lI) <= 0.f) {
                zero_weight_layer=true;
                pathR = pathI = 0.f;
                break;
            }
            current_state = next_state;

        }
        if (!skip_path || !is_zero || !zero_weight_layer) {
            sumR += pathR;
            sumI += pathI;
        }

        path_counter++;
        if (!zero_weight_layer) path_NZ_counter++;

        // compute next path
        // updating ndxs[]
        int ll;
        for (ll=(((is_zero || zero_weight_layer) && l<(L-1))? l : L-2); ll>=0 ; ll--) {
            ndxs[ll]++;
            start_layer=ll;
            if (ndxs[ll]==N && ll>0)  { 
                ndxs[ll] = 0;
            }
            else
                break;
        }
    } // main simulation loop (while)
    aR = sumR;
    aI = sumI;

    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu paths, %llu non zero\n", path_counter, path_NZ_counter);
}
