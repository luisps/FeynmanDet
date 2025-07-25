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

#include "simulator_all.hpp"

#if defined(_OPENMP)
#include <omp.h>
extern int n_threads;
#endif

void simulate_all_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    const int NQ=circuit->size->num_qubits;
    const StateT N = 1 << NQ;
    StateT path_counter=0, path_NZ_counter=0;

    double const total_paths = pow(2.F, (double)(NQ*(L-1)));
    printf ("%le existing paths\n", total_paths);
    
    /*int init_state_arr[NQ];
    int final_state_arr[NQ];
    for (int i=0; i<NQ; i++){
        init_state_arr[i]=qb_value(i,init_state);
        final_state_arr[i]=qb_value(i,final_state);	
    }*/
    aR = aI = 0.f;

    // Simulation starts

    // omp parallel block
#if defined(_OPENMP)
    if (n_threads==-1) {
        int const n_cores = omp_get_num_procs();
        n_threads = (N>n_cores ? n_cores : N);
    }
    fprintf (stdout, "OpenMP: %d threads\n", n_threads);
#pragma omp parallel num_threads(n_threads) proc_bind(spread)
#endif
    {
        StateT path_counterL=0, path_NZ_counterL=0;
        // main loop
        StateT ndxs0;
        float sumR=0.f, sumI=0.f;
#if defined(_OPENMP)
        double T_totaltime=0.;
#pragma omp for schedule(dynamic,1)
#endif
        for (ndxs0 = 0 ; ndxs0 < N ; ndxs0++) {
#if defined(_OPENMP)
            clock_t start, end;
            start=clock();
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
            for (int i=1 ; i<L-1 ; i++) ndxs[i]=0 ;
            
            float wR[L-1], wI[L-1];
            
            // early termination if the amplitude
            // from init_state to ndxs[0] is zero
            float pathR = 1.f;
            float pathI = 0.f;
            TCircuitLayer *layer = &circuit->layers[0];
            layer_w(layer, 0, init_state, ndxs0, pathR, pathI);
            
            
            if (complex_abs_square(pathR, pathI) <= 0.f) {
                continue;
            }
            wR[0]=pathR;
            wI[0]=pathI;
            int start_layer=1;

            while (ndxs[1] < N) {
                
                /*
                 fprintf (stderr, "ndxs= ");
                 fprintf (stderr, "%llu ", ndxs0);
                 for (int lll=1; lll<L-1 ; lll++)
                 fprintf (stderr, "%llu ", ndxs[lll]);
                 fprintf (stderr, "\n");*/
                
                pathR = (start_layer==0? 1.f : wR[start_layer-1]);
                pathI = (start_layer==0? 0.f : wI[start_layer-1]);
                StateT current_state = (start_layer==0? init_state : ndxs[start_layer-1]);
                
                int l;
                StateT next_state;
                bool zero_weight_layer=false;
                
                // iterate over layers
                for (l=start_layer ; l<L ; l++) {
                    float lR=1.f;
                    float lI=0.f;
                    next_state = (l< L-1 ? ndxs[l] : final_state);
                    
                    layer = &circuit->layers[l];
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
                
                if (!zero_weight_layer) {
                    sumR += pathR;
                    sumI += pathI;
                    
                }
                path_counterL++;
                /*if (!(path_counter & 0x00FFFFF)) {
                 fprintf(stderr, "\rpath_counter=%llu", path_counter);
                 }*/
                if (!zero_weight_layer) path_NZ_counterL++;
                
                // compute next path
                // if l==0 (0 amplitude in the first layer)
                // break off from the while loop
                // (to iterate over ndxs0)
                if (l==0) break;
                // updating ndxs[1..L-2]
                int ll;
                for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>0 ; ll--) {
                    ndxs[ll]++;
                    start_layer=ll;
                    if (ndxs[ll]==N && ll>1)  {
                        ndxs[ll] = 0;
                    }
                    else
                        break;
                }
            } // main simulation loop (while)
#if defined(_OPENMP)
            end=clock();
            double time_taken=double(end - start)/double(CLOCKS_PER_SEC);
            fprintf (stdout, "thread %d with ndxs0=%llu took %lf secs\n", omp_get_thread_num(), ndxs0, time_taken);
            T_totaltime += time_taken;
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
#if defined(_OPENMP)
        printf ("Thread %d: %llu evaluated paths, %llu non zero (%lf secs)\n", omp_get_thread_num(), path_counterL, path_NZ_counterL, T_totaltime);

#endif
    } // end omp parallel

    printf ("\n");
    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu evaluated paths, %llu non zero\n", path_counter, path_NZ_counter);
}
