//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#include <math.h>
#include <algorithm>
#include "my_complex.h"
#include "layer.hpp"
#include <time.h>
#include <iostream>

#include "colouring_RG.hpp"

#include "simulator_PB.hpp"

#if defined(_OPENMP)
#include <omp.h>
extern int n_threads;
#endif

static void print_PB (colourT *states, int NQ, int L) {
    fprintf (stderr, "Colouring: \n");
    for (int q = 0; q < NQ; q++) {
        fprintf (stderr, "qb%d: ", q);
        for (int l = 1; l < L; l++) {
            colourT colour= states[q + (l - 1)*NQ];
            char s[5];
            switch (colour & 0x0FF) {
                case INVALID:
                    snprintf (s,4,"I");
                    break;
                case GREEN0:
                    snprintf (s,4,"G0");
                    break;
                case GREEN1:
                    snprintf (s,4,"G1");
                    break;
                case RED:
                    snprintf (s,4,"R");
                    break;
                case PINK:
                    snprintf (s,4,"P");
                    break;
                case BLUE:
                    snprintf (s,4,"B");
                    break;
                case BLUE_I:
                    snprintf (s,4,"BI");
                    break;
                case BLUE_X:
                    snprintf (s,4,"BX");
                    break;
                case BLUE_CX:
                    snprintf (s,4,"BCX");
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



// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
// gates 2 qubits
// G2P0 - 2 qubits, no parameters
// 'id2'  - 20    Used for errors
// 'cx'  - 21 B
// 'cz'  - 22
// G2P1 - 2 qubits, 1 parameter
// 'cp'  - 31



void printBits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    
    for (i = size-1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

// verify whether this state is allowed by colouring
bool validate_PB (StateT const state, StateT const prev_state, int const *const colours, int const NQ) {
    bool valid = true;

    
    for (int i=0; i<NQ && valid ; i++){
        int const next_state_CL = colours[i];
        
        if (next_state_CL== GREEN0 || next_state_CL== GREEN1) {
            if (next_state_CL != qb_value(i,state)) {
                valid = false; // break from inner loop
            }
        }
        else if (next_state_CL == BLUE_I || next_state_CL == BLUE_X) {
            int const pred_qb_value = (next_state_CL == BLUE_I ? qb_value(i,prev_state) : !qb_value(i,prev_state));
            if (pred_qb_value != qb_value(i,state)) {
                valid=false;
            }
        }
        /*else if (next_state_CL == BLUE) {
            int const pred_qb_value = qbv_given_previous_gate(circuit, ll, i, prev_state);
            if (pred_qb_value != qb_value(i,ndxs[ll])) {
                invalid_state_green=true;
            }
        }*/
        else if (next_state_CL > 255 || next_state_CL == BLUE_CX) {
            // get the control qubit index
            int const c_qb = next_state_CL >> 8;
            // get the control qubit value
            int const c_qb_v = qb_value(c_qb, prev_state);
            int const pred_qb_value = (c_qb_v == 0 ? qb_value(i,prev_state) : !qb_value(i,prev_state));
            if (pred_qb_value != qb_value(i,state)) {
                valid = false;
            }
        }

    } // for to verify qubits and green
    return valid;
}

void simulate_PB_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    const int NQ=circuit->size->num_qubits;
    const StateT N = 1 << NQ;
    StateT path_counter=0, path_NZ_counter=0;

    if (L<4) { // 4 layers are required
        fprintf (stderr, "The circuit has %d layers: 4 is the minimum!\n", L);
        return ;
    }
    
    double const total_paths = pow(2.F, (double)(NQ*(L-1)));
    fprintf (stdout, "%le existing paths\n", total_paths);
    
    fflush(stderr);
    fflush(stdout);



    int init_state_arr[NQ];
    int final_state_arr[NQ];
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

    aR = aI = 0.f;
    
    // verify if there is a -1 in the colouring
    // if yes, amplitude = 0 + 0j
    for(int s=0; s<(L-1)*NQ; s++){
        if(colours[s]==INVALID){
            printf ("EARLY TERMINATION: < %llu | U | %llu > = %.6f + i %.6f\n", (unsigned long long)final_state,
                (unsigned long long)init_state, aR, aI);
            return;
        }
    }

    fs_PB(circuit, colours);
    
    //print_PB (colours, NQ, L);
    
    // Simulation starts
    // omp parallel block
#if defined(_OPENMP)
    if (n_threads==-1) {
        int const n_cores = omp_get_num_procs();
        n_threads = (N>n_cores ? n_cores : N);
    }
    fprintf (stderr, "OpenMP: %d threads required\n", n_threads);
#pragma omp parallel num_threads(n_threads) proc_bind(spread)
    {
        int const threadID = omp_get_thread_num();
        
#pragma omp single nowait
        {
        fprintf (stderr, "OpenMP: thread %d reports %d threads\n", threadID, omp_get_num_threads());
        fflush (stderr);
        }
#else
    {
        int const threadID = -1;
        fprintf (stderr, "OpenMP NOT IDENTIFIED\n");
        fflush (stderr);
#endif

        StateT path_counterL=0, path_NZ_counterL=0;
    
        // explicitly include up to 2 for loops
        // for parallel execution by OpenMP
        
        // We will be using openMP parallel for
        // and eventually collapsing 2 (two) for loops
        // The next pre-processor constant tells us whetehr or not
        // to collapse
        StateT ndxs0;

        float sumR=0.f, sumI=0.f;
#if defined(_OPENMP)
//#define _COLLAPSE2
#if defined(_COLLAPSE2)
        StateT ndxs1;
#endif
        double start, end;
        start=omp_get_wtime();
        int n_tasks=0;
#if defined(_COLLAPSE2)
#pragma omp for schedule(dynamic) collapse(2)
        for (ndxs0 = 0 ; ndxs0 < N ; ndxs0++) {
            for (ndxs1 = 0 ; ndxs1 < N ; ndxs1++) {
#else
#pragma omp for schedule(dynamic)
                for (ndxs0 = 0 ; ndxs0 < N ; ndxs0++) {
#endif
#else
                    for (ndxs0 = 0 ; ndxs0 < N ; ndxs0++) {
#endif
                        
#if defined(_OPENMP)
                        n_tasks++;
#endif
                        
                        
                        // all intermediate layers indexes to 0
                        // except intermediate layer 0
                        // this one iterates as a for loop, to facilitate OpenMP
                        StateT ndxs[L-1];
                        // we make ndxs[0] equal to ndxs0
                        // only to avoid conditionals below
                        // this is only for reading
                        // all writes must be made to ndxs0
                        ndxs[0] = ndxs0;
#if defined(_COLLAPSE2)
                        ndxs[1] = ndxs1;
#endif
#if defined(_COLLAPSE2)
                        int const Collapsed_loops=2;
#else
                        int const Collapsed_loops=1;
#endif
                        for (int i=Collapsed_loops ; i<L-1 ; i++) ndxs[i]=0 ;
                        
                        float wR[L-1], wI[L-1];
                        
                        float pathR = 1.f;
                        float pathI = 0.f;
                        float lR=1.f;
                        float lI=0.f;
                        
                        // We have to validate whether ndxs0 is a valid state
                        // given the colouring
                        if(!validate_PB (ndxs0, init_state, &colours[0*NQ], NQ)){
                            
                            continue;  // if not valid and _COLLAPSE2 skip to next ndxs1
                            // very unfortunate
                        }
                        
                        // early termination if the amplitude
                        // from init_state to ndxs[0] is zero
                        TCircuitLayer *layer = &circuit->layers[0];
                        layer_w(layer, 0, init_state, ndxs0, lR, lI);
                        
                        // Get this very carefully
                        // only proceed with this iteration
                        // which in fact is a pair (ndxs0,ndxs1)
                        // if the amplitude wasn't zero
                        if (complex_abs_square(lR, lI) <= 0.f) {
                            continue;   // if zero and _COLLAPSE2 skip to next ndxs1
                            // very unfortunate
                        }
                        
                        wR[0]=pathR = lR;
                        wI[0]=pathI = lI;
                        
#if defined(_COLLAPSE2)
                        // We have to validate whether ndxs1 is a valid state
                        // given the colouring
                        if(!validate_PB (ndxs1, ndxs0, &colours[1*NQ], NQ)){
                            
                            continue;  // if not valid and _COLLAPSE2 skip to next ndxs1
                        }
                        // early termination if the amplitude
                        // from ndxs[0] to ndxs[1] is zero
                        lR=1.f; lI=0.f;
                        layer = &circuit->layers[1];
                        layer_w(layer, 1, ndxs0, ndxs1, lR, lI);
                        if (complex_abs_square(lR, lI) <= 0.f) {
                            continue;  // next ndxs1
                        }
                        complex_multiply(pathR, pathI, lR, lI, pathR, pathI);
                        wR[1]=pathR;
                        wI[1]=pathI;
#endif
                        
                        int start_layer=Collapsed_loops;
                        
                        // main simulation loop
                        while (ndxs[Collapsed_loops] < N) {
                            
                            pathR = (start_layer==0? 1.f : wR[start_layer-1]);
                            pathI = (start_layer==0? 0.f : wI[start_layer-1]);
                            StateT current_state = (start_layer==0? init_state : ndxs[start_layer-1]);
                            
                            int l;
                            StateT next_state;
                            bool zero_weight_layer=false;
                            
                            // iterate over layers
                            for (l=start_layer ; l<L ; l++) {
                                lR=1.f; lI=0.f;
                                next_state = (l< L-1 ? ndxs[l] : final_state);
                                
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
                            if  (!zero_weight_layer) {
                                sumR += pathR;
                                sumI += pathI;
                                path_NZ_counterL++;
                                
                                // DEBUG
                                /*printf ("Non zero path: ");
                                 for (int lll=0 ; lll<L-1 ; lll++) {
                                 printf ("%llu ", ndxs[lll]);
                                 }
                                 printf ("= %e + i %e\n", pathR, pathI);*/
                                
                            }
                            path_counterL++;
                            
                            // compute next path
                            // if l==0 (0 amplitude in the first layer)
                            // break off from the while loop
                            // (to iterate over ndxs1)
                            if (l==Collapsed_loops-1) break;
                            
                            // updating ndxs[Collapsed_loops..L-2]
                            // compute next path skipping invalid GREENs
                            
                            // compute next path skipping invalid GREENs and fixed BLUEs
                            int ll;
                            bool invalid_state_green = true;
                            for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=Collapsed_loops && invalid_state_green ; ll--) {
                                
                                //printf ("change ndxs[%d]=%llu\n", ll, ndxs[ll]);
                                // get the state from the previous layer. Will need it
                                StateT const prev_state = (ll==0 ? init_state : ndxs[ll-1]);
                                
                                invalid_state_green = true;
                                while(invalid_state_green) {
                                    ndxs[ll]++;
                                    //if (ll==0)  fprintf (stderr, "ndxs[0] = %llu \n", ndxs[0]);
                                    
                                    if (ndxs[ll]==N && ll>Collapsed_loops)  { // this layer overflows
                                        ndxs[ll] = 0;
                                        break;        // break only from inner loop
                                    }
                                    else if (ndxs[Collapsed_loops]==N && ll==Collapsed_loops)  { // back to for loops
                                        invalid_state_green = false; // terminate outer loop
                                        break;   // terminate inner loop
                                    }
                                    start_layer=ll;
                                    //printf ("changed ndxs[%d]=%llu (ndxs[0] = %llu) \n", ll, ndxs[ll], ndxs[0]);
                                    
                                    // verify whether this ndxs complies with the colouring
                                    invalid_state_green = !validate_PB(ndxs[ll], prev_state, &colours[ll*NQ], NQ);
                                    /*if (invalid_state_green) { // need to know invalid bit index
                                     ndxs[ll] |= ((1 << i)- 1) ; // skip all  intermediate non valid states
                                     }*/
                                    
                                }  // while (invalid_state_green)
                                //printf ("END FOR LOOP ndxs[%d]=%llu\n", ll, ndxs[ll]);
                            }    // for backward change layers ndxs
                            
                        } // main simulation loop (while)
#if defined(_COLLAPSE2)
                    }  // main simulation loop (ndxs1 and omp for)
#endif
                }  // main simulation loop (ndxs0 and omp for)
#if defined(_OPENMP)
#pragma omp atomic
#endif
                aR += sumR;
#if defined(_OPENMP)
#pragma omp atomic
#endif
                aI += sumI;
#if defined(_OPENMP)
#pragma omp atomic
#endif
                path_counter += path_counterL;
#if defined(_OPENMP)
#pragma omp atomic
#endif
                path_NZ_counter += path_NZ_counterL;
                /*
                 #if defined(_OPENMP)
                 end=omp_get_wtime();
                 double time_taken=double(end - start)*1000.F;
                 printf ("Thread %d: %llu evaluated paths, %llu non zero (%.2lf mili secs), n_tasks=%d\n", omp_get_thread_num(), path_counterL, path_NZ_counterL, time_taken, n_tasks);
                 
                 #endif
                 */
            }// end omp parallel
            
            printf ("\n");
    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu paths, %llu non zero\n", path_counter, path_NZ_counter);
}
