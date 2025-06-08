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
#endif

void simulate_all_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    StateT *ndxs = new StateT[L-1];
    const int NQ=circuit->size->num_qubits;
    const StateT N = 1 << NQ;
    StateT path_counter=0, path_NZ_counter=0;

    float* wR=new float[L-1];
    float* wI=new float[L-1];

    double const total_paths = pow(2.F, (double)(NQ*(L-1)));
    printf ("%le existing paths\n", total_paths);
    
    int init_state_arr[NQ];
    int final_state_arr[NQ];
    float sumR=0.f, sumI=0.f;
    for (int i=0; i<NQ; i++){
        init_state_arr[i]=qb_value(i,init_state);
        final_state_arr[i]=qb_value(i,final_state);	
    }

    int start_layer=0;

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
        
        float pathR = (start_layer==0? 1.f : wR[start_layer-1]);
        float pathI = (start_layer==0? 0.f : wI[start_layer-1]);
        StateT current_state = (start_layer==0? init_state : ndxs[start_layer-1]);

        int l;
        StateT next_state;
        bool zero_weight_layer=false;

        // iterate over layers
        for (l=start_layer ; l<L ; l++) {
            float lR=1.f;
            float lI=0.f;
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
        
        }
            
        if (!zero_weight_layer) {
            sumR += pathR;
            sumI += pathI;
        }
        path_counter++;
        if (!(path_counter & 0x00FFFFF)) {
            fprintf(stderr, "\rpath_counter=%llu", path_counter);
        }
        if (!zero_weight_layer) path_NZ_counter++;
        
        // compute next path
        // updating ndxs[]
        int ll;
        for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=0 ; ll--) {
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

    printf ("\n");
    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu evaluated paths, %llu non zero\n", path_counter, path_NZ_counter);
}
