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

//#define _OPENMP

#if defined(_OPENMP)
#include <omp.h>
extern int n_threads;

#define CHUNKSIZE 1

#endif

void simulate_all_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    const int NQ=circuit->size->num_qubits;
    const StateT N = 1 << NQ;
    StateT path_counter=0, path_NZ_counter=0;
    
    if (L<4) { // 4 layers are required
        fprintf (stderr, "The circuit has %d layers: 4 is the minimum!\n", L);
        fflush(stderr);
        return ;
    }

    double const total_paths = pow(2.F, (double)(NQ*(L-1)));
    fprintf (stdout, "%le existing paths\n", total_paths);
    
    fflush(stdout);
    
    aR = aI = 0.f;

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

        // explicitly include up to 2 for loops
        // for parallel execution by OpenMP
        
        // We will be using openMP parallel for
        // and eventually collapsing 3 (three) for loops
        // The next pre-processor constant tells us whetehr or not
        // to collapse
        StateT ndxs0;

#if defined(_OPENMP)
#define _SCRAMBLE
#if defined(_SCRAMBLE)
        const uint64_t maskN = N - 1;
        const int m = __builtin_ctzll(N);   // since N is power of two

        // Choose affine permutation parameters
        const uint64_t a = 6364136223846793005ULL;  // must be odd
        const uint64_t b = 1442695040888963407ULL;

#endif
        double start, end;
        start=omp_get_wtime();
        int n_tasks=0;

#if defined(_SCRAMBLE)
#pragma omp for schedule(static, CHUNKSIZE)
        for (StateT t = 0 ; t < N ; t++) {
            ndxs0 = (StateT) ((a * t + b) & maskN);
#else  // NO_SCRAMBLE
#pragma omp for schedule(static, CHUNKSIZE)
        for (ndxs0 = 0 ; ndxs0 < N ; ndxs0++) {
#endif   // end _SCRAMBLE
#else   // NO _OPENMP
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

            for (int i=1 ; i<L-1 ; i++) ndxs[i]=0 ;
                                
            float wR[L-1], wI[L-1];
                                
            float pathR = 1.f;
            float pathI = 0.f;
            float lR=1.f;
            float lI=0.f;
                                
            // early termination if the amplitude
            // from init_state to ndxs[0] is zero
            TCircuitLayer *layer = &circuit->layers[0];
            layer_w(layer, 0, init_state, ndxs0, lR, lI);
                                
            // only proceed with this iteration
            // if the amplitude wasn't zero
            if (complex_abs_square(lR, lI) <= 0.f) {
                continue;  // next t
            }
                                
            wR[0]=pathR = lR;
            wI[0]=pathI = lI;
                                
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
                    lR=1.f; lI=0.f;
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
                } // end iterating layers
                                    
                if (!zero_weight_layer) {
                    sumR += pathR;
                    sumI += pathI;
                    path_NZ_counterL++;
                                        
                    //fprintf (stderr, "(T%d) - %llu th NZ path (%llu -> %llu -> %llu ...\n", threadID, path_NZ_counterL, init_state, ndxs0, ndxs1);
                    //fflush (stderr);
                }
                path_counterL++;
                /*if (!(path_counter & 0x00FFFFF)) {
                    fprintf(stderr, "\rpath_counter=%llu", path_counter);
                }*/
                                    
                // compute next path
                // if l==0 (0 amplitude in the first layer)
                // break off from the while loop
                // (to iterate over t)
                if (l==0) break;
                // updating ndxs[1..L-2]
                int ll;
                for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=1 ; ll--) {
                    ndxs[ll]++;
                    start_layer=ll;
                    if (ndxs[ll]==N && ll>1)  {
                        ndxs[ll] = 0;
                    }
                    else
                        break;
                }
            } // main simulation loop (while)
        }  // main simulation loop (t and omp for)
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
        end=omp_get_wtime();
        double time_taken=double(end - start)*1000.F;
        printf ("Thread %d: %llu evaluated paths, %llu non zero (%.2lf mili secs), n_tasks=%d\n", omp_get_thread_num(), path_counterL, path_NZ_counterL, time_taken, n_tasks);

#endif
         
    } // end omp parallel

    printf ("\n");
    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu evaluated paths, %llu non zero\n", path_counter, path_NZ_counter);
}
