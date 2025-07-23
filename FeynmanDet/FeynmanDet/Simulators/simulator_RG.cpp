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

static void print_RG (colourT *states, int NQ, int L) {
    fprintf (stderr, "Colouring: \n");
    for (int q = 0; q < NQ; q++) {
        fprintf (stderr, "qb%d: ", q);
        for (int l = 1; l < L; l++) {
            colourT colour= states[q + (l - 1)*NQ];
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
                default:
                    snprintf (s,4,"ERR");
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
    colourT* colours;
    colours=fs_bs_RG(circuit, init_state_arr, final_state_arr);

    print_RG (colours, NQ, L);
    int start_layer=0;

    // verify if there is a -1 in the colouring
    // if yes, amplitude = 0 + 0j
    for(int s=0; s<(L-1)*NQ; s++){
        if(colours[s]==INVALID){
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

        //fprintf (stderr, "Main LOOP iteration\n");
        /*
        fprintf (stderr, "ndxs= ");
        for (int lll=0; lll<L-1 ; lll++)
            fprintf (stderr, "%llu ", ndxs[lll]);
        fprintf (stderr, "\n");*/

        float pathR = (start_layer==0? 1.f : wR[start_layer-1]);
        float pathI = (start_layer==0? 0.f : wI[start_layer-1]);
        StateT current_state = (start_layer==0? init_state : ndxs[start_layer-1]);

        int l;
        StateT next_state;
        bool zero_weight_layer=false;
        //bool is_zero=false;
        
        // iterate over layers
        for (l=start_layer ; l<L ; l++) {
            float lR=1.f;
            float lI=0.f;
            next_state = (l< L-1 ? ndxs[l] : final_state);
            /*is_zero = false;

            // verify whether this next state is allowed by colouring
            if(l < L-1){
                for (int i=0; i<NQ; i++){
                    
                    int const next_state_CL = colours[i+l*NQ];
                    if (next_state_CL==2) continue;
                    
                    if (next_state_CL != qb_value(i,next_state)) {
                        pathR = pathI = 0.f;
                        is_zero=true;
                        break;
                    }
                }
            }
            if(is_zero){
                break;
            }*/
        
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

        } // end iterating layers
        //if  (!is_zero && !zero_weight_layer) {
        if  (!zero_weight_layer) {
            sumR += pathR;
            sumI += pathI;
            path_NZ_counter++;
            
            // DEBUG
            /*printf ("Non zero path: ");
            for (int lll=0 ; lll<L-1 ; lll++) {
                printf ("%llu ", ndxs[lll]);
            }
            printf ("= %e + i %e\n", pathR, pathI);*/

        }
        path_counter++;

        // compute next path
        // updating ndxs[]
        
        // all paths
        /*int ll;
        //for (ll=(((is_zero || zero_weight_layer) && l<(L-1))? l : L-2); ll>=0 ; ll--) {
        for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=0 ; ll--) {
            ndxs[ll]++;
            start_layer=ll;
            if (ndxs[ll]==N && ll>0)  { 
                ndxs[ll] = 0;
            }
            else
                break;
        }*/
        
        
        // compute next path skipping invalid GREENs
        int ll;
        bool invalid_state_green = true;
        for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=0 && invalid_state_green ; ll--) {

            //printf ("change ndxs[%d]=%llu\n", ll, ndxs[ll]);
            invalid_state_green = true;
            while(invalid_state_green) {
                ndxs[ll]++;
                if (ndxs[ll]==N && ll>0)  { // this layer overflows
                    ndxs[ll] = 0;
                    break;        // break only from inner loop
                }
                else if (ndxs[0]==N && ll==0)  { // simulation finished
                    invalid_state_green = false; // terminate outer loop
                    break;   // terminate inner loop
                }
                start_layer=ll;
                //printf ("changed ndxs[%d]=%llu\n", ll, ndxs[ll]);
                invalid_state_green = false;
                // verify whether this ndxs complies with the colouring
                for (int i=0; i<NQ && !invalid_state_green ; i++){
                    
                    int const next_state_CL = colours[i+ll*NQ];
                    
                    if (next_state_CL== GREEN0 || next_state_CL== GREEN1) {
                        if (next_state_CL != qb_value(i,ndxs[ll])) {
                            invalid_state_green=true; // break from all loops
                            //ndxs[ll] += (1 << i)- 1 ; // skip all  intermediate non valid states
                        }
                    }
                } // for to verify qubits and green
            }  // while (invalid_state_green)
            //printf ("END FOR LOOP ndxs[%d]=%llu\n", ll, ndxs[ll]);
        }    // for backward change layers ndxs

    } // main simulation loop (while)
    fprintf (stderr, "Main LOOP terminated\n");
    aR = sumR;
    aI = sumI;

    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu paths, %llu non zero\n", path_counter, path_NZ_counter);
}
