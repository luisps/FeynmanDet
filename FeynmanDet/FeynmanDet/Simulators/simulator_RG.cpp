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

#if defined(_OPENMP)
#include <omp.h>
extern int n_threads;
#endif

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

// verify whether this state is allowed by colouring
// return:  the index of the violating qubit
//          -1 if valid
int validate_RG (StateT const state, int const *const colours, int const NQ) {
    int invalid = -1;

    for (int i=0; i<NQ && invalid==-1; i++){
        
        int const next_state_CL = colours[i];
        
        if (next_state_CL==2) continue;  // RED is OK, accepts all
        else if (next_state_CL != qb_value(i,state)) { // G0!= 1 or G1!= 0
            invalid=i;
        }
    }
    return invalid;
}

void simulate_RG_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

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


    // Colouring
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

    //print_RG (colours, NQ, L);

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
            
            float sumR=0.f, sumI=0.f;
#if defined(_OPENMP)
            // explicitly include up to 2 for loops
            // for parallel execution by OpenMP
            
            // We will be using openMP parallel for
            // and eventually collapsing 2 (two) for loops
            // The next pre-processor constant tells us whetehr or not
            // to collapse
            StateT ndxs0;
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
#if defined(_OPENMP)
                        ndxs[0] = ndxs0;
                        int const Collapsed_loops=1;
#if defined(_COLLAPSE2)
                        ndxs[1] = ndxs1;
                        int const Collapsed_loops=2;
#endif
#else
                        int const Collapsed_loops=0;
#endif
                        for (int i=Collapsed_loops ; i<L-1 ; i++) ndxs[i]=0 ;
                        
                        float wR[L-1], wI[L-1];
                        
                        float pathR = 1.f;
                        float pathI = 0.f;
                        float lR=1.f;
                        float lI=0.f;
                        
#if defined(_OPENMP)
                        // We have to validate whether ndxs0 is a valid state
                        // given the colouring
                        if(validate_RG (ndxs0, &colours[0*NQ], NQ)!=-1){
                            
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
                        if(validate_RG (ndxs1, &colours[1*NQ], NQ)!=-1){
                            
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
#endif
                        
                        int start_layer=Collapsed_loops;
                        
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
                            
                            int ll;
                            int invalid_state_qb = 0;  // -1 is valid
                            for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=Collapsed_loops && invalid_state_qb != -1 ; ll--) {
                                
                                invalid_state_qb = 0;  // -1 is valid
                                while(invalid_state_qb!=-1) {
                                    ndxs[ll]++;
                                    if (ndxs[ll]==N && ll>Collapsed_loops)  { // this layer overflows
                                        ndxs[ll] = 0;
                                        break;        // break only from inner loop
                                    }
                                    else if (ndxs[Collapsed_loops]==N && ll==Collapsed_loops)  { // back to for loops
                                        invalid_state_qb = -1; // terminate outer loop
                                        break;   // terminate inner loop
                                    }
                                    start_layer=ll;
                                    //printf ("changed ndxs[%d]=%llu\n", ll, ndxs[ll]);
                                    // verify whether this ndxs complies with the colouring
                                    invalid_state_qb = validate_RG(ndxs[ll], &colours[ll*NQ], NQ);
                                    if (invalid_state_qb != -1) { // need to know invalid bit index
                                        ndxs[ll] |= ((1 << invalid_state_qb)- 1) ; // skip all  intermediate non valid states
                                    }
                                }  // while (invalid_state_green)
                                //printf ("END FOR LOOP ndxs[%d]=%llu\n", ll, ndxs[ll]);
                            }  // for backward change layers ndxs
                            
                        } // main simulation loop (while)
#if defined(_OPENMP)
#if defined(_COLLAPSE2)
                    }  // main simulation loop (ndxs1 and omp for)
#endif
                }  // main simulation loop (ndxs0 and omp for)
#endif
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
            }// end omp parallel (NOTE: { included even if !NOT _OPENMP)
            printf ("\n");
            printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
            printf ("%llu paths, %llu non zero\n", path_counter, path_NZ_counter);
        }  // end function
            
            void simulate_RG_paths_old (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

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


                // Colouring
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

                //print_RG (colours, NQ, L);

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
                                    if(validate_RG (ndxs0, &colours[0*NQ], NQ)!=-1){
                                        
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
                                    if(validate_RG (ndxs1, &colours[1*NQ], NQ)!=-1){
                                        
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
                                        
                                        int ll;
                                        int invalid_state_qb = 0;  // -1 is valid
                                        for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=Collapsed_loops && invalid_state_qb != -1 ; ll--) {
                                            
                                            invalid_state_qb = 0;  // -1 is valid
                                            while(invalid_state_qb!=-1) {
                                                ndxs[ll]++;
                                                if (ndxs[ll]==N && ll>Collapsed_loops)  { // this layer overflows
                                                    ndxs[ll] = 0;
                                                    break;        // break only from inner loop
                                                }
                                                else if (ndxs[Collapsed_loops]==N && ll==Collapsed_loops)  { // back to for loops
                                                    invalid_state_qb = -1; // terminate outer loop
                                                    break;   // terminate inner loop
                                                }
                                                start_layer=ll;
                                                //printf ("changed ndxs[%d]=%llu\n", ll, ndxs[ll]);
                                                // verify whether this ndxs complies with the colouring
                                                invalid_state_qb = validate_RG(ndxs[ll], &colours[ll*NQ], NQ);
                                                if (invalid_state_qb != -1) { // need to know invalid bit index
                                                 ndxs[ll] |= ((1 << invalid_state_qb)- 1) ; // skip all  intermediate non valid states
                                                 }
                                            }  // while (invalid_state_green)
                                            //printf ("END FOR LOOP ndxs[%d]=%llu\n", ll, ndxs[ll]);
                                        }  // for backward change layers ndxs
                                        
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

