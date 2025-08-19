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

#include "simulator_RG_new_errors.hpp"

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

void simulate_RG_paths_new_errors (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    // current state at each layer
    StateT *state = new StateT[L-1];
    // state iterators
    FixedBitsSequence2 *state_it = new FixedBitsSequence2[L-1];
    const int NQ=circuit->size->num_qubits;
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


    // init the states and verify if there is a -1 in the colouring
    // if yes, amplitude = 0 + 0j
    for (int l = 0; l < L-1; l++) {
        int n_fixed_bits;
        int fixed_bits[sizeof(StateT)];
        int fixed_values[sizeof(StateT)];
        n_fixed_bits=0;
        for (int q = 0; q < NQ; q++) {
            colourT colour= colours[q + l*NQ];
            switch (colour) {
                case INVALID:
                    aR=0.f;
                    aI=0.f;
                    printf ("EARLY TERMINATION: < %llu | U | %llu > = %.6f + i %.6f\n", (unsigned long long)final_state,
                            (unsigned long long)init_state, aR, aI);
                    return;
                    break;
                case GREEN0:
                case GREEN1:
                    fixed_bits[n_fixed_bits] = q;
                    fixed_values[n_fixed_bits] = colour;
                    n_fixed_bits++;
                    break;
            }
        }  // end iterate over this layer qubits
        // initialize iterator
        state_it[l].init_iterator(NQ, n_fixed_bits, fixed_bits, fixed_values);
    }

    // Simulation starts

    // all intermediate layers are placed in first in sequence
    for (int l=0 ; l<L-1 ; l++) {
        state_it[l].generate_next_in_sequence(state[l]) ;
    }
    
    // main simulation loop
    bool finish_simulation=false;
    while ( ! finish_simulation)  {

        //fprintf (stderr, "Main LOOP iteration\n");
        /*
        fprintf (stderr, "ndxs= ");
        for (int lll=0; lll<L-1 ; lll++)
            fprintf (stderr, "%llu ", ndxs[lll]);
        fprintf (stderr, "\n");*/

        float pathR = (start_layer==0? 1.f : wR[start_layer-1]);
        float pathI = (start_layer==0? 0.f : wI[start_layer-1]);
        StateT current_state = (start_layer==0? init_state : state[start_layer-1]);

        int l;
        StateT next_state;
        bool zero_weight_layer=false;
        //bool is_zero=false;
        
        // iterate over layers
        for (l=start_layer ; l<L ; l++) {
            float lR=1.f;
            float lI=0.f;
            next_state = (l< L-1 ? state[l] : final_state);
        
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
                printf ("%llu ", state[lll]);
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
        
        for (int ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=0 ; ll--) {
            
            if (ll >0 && !state_it[ll].generate_next_in_sequence(state[ll]))  { // this layer overflows, reset and go to previous layer
                state_it[ll].reset();
                state_it[ll].generate_next_in_sequence(state[ll]) ;
            }
            else if (ll==0 && !state_it[0].generate_next_in_sequence(state[0]))  { // simulation finished
                finish_simulation = true;
                break;   // terminate iterating states
            }
            else {  // continue on this layer ; break out of for loop
                break;
            }
        }

    } // main simulation loop (while)
    fprintf (stderr, "Main LOOP terminated\n");
    aR = sumR;
    aI = sumI;

    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu paths, %llu non zero\n", path_counter, path_NZ_counter);
}
