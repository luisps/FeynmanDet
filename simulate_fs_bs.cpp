//
//  simulate_all.cpp
//  Path_simulator
//
//  Created by Luis Paulo Santos on 03/11/2023.
//
#include "fs_bs_working.cpp"
#include "simulate_all.hpp"
#include <math.h>
#include "complex.h"
#include "layer.hpp"


void simulate_all_paths (TCircuit *circuit,
                         int init_state, int final_state,
                         float& aR, float& aI) {
    
    const int L = circuit->size->num_layers;
    int *ndxs = new int[L-1];
    const int N = (int)(powf(2.f, (float)circuit->size->num_qubits));

    
    float sumR=0.f, sumI=0.f;
    
    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;
    
    while (ndxs[0] < N) {

        // .... process path : init_state-> ndxs[0] -> ..... -> ndxs[L-2] -> final_state
        
        float pathR=1.f, pathI=0.f;
                
        int current_state = init_state;
        int next_state, l;
        
        for (l=0 ; l<L ; l++) {
            float lR=1.f, lI=0.f;
            
            next_state = (l< L-1 ? ndxs[l] : final_state);
            
            // processar o layer (lR, lI)
            // get gate layer l
            TCircuitLayer *layer = &circuit->layers[l];
            
            // get the amplitude of this layer for the given current and next states
            // Wer are reevaluating all layers even though a single ndxs[] changes: INNEFICIENT
            layer_w (layer, l, current_state, next_state, lR, lI);
            
            if (complex_abs_square(lR, lI) <= 0.) {
                ///// TERMINATE PATH
                pathR = pathI = 0.f;
                continue;
            }

            complex_multiply (pathR, pathI, lR, lI, pathR, pathI);
            
            current_state = next_state;

        }
        
        sumR += pathR;
        sumI += pathI;

        
        // iterate to next path
        int ll;
        for (ll=L-2 ; ll>=0 ; ll--) {
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
    
}




void simulate_all_paths_opt (TCircuit *circuit,
                         int init_state, int final_state,
                         float& aR, float& aI) {
    
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits; 
    int *ndxs = new int[L-1];
    const int N = (int)(powf(2.f, (float)NQ));
    
    float sumR=0.f, sumI=0.f;
    
    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;
    
    // store the amplitudes per layer
    float *layer_aR = new float[L];
    float *layer_aI = new float[L];

    // compute the initial layers' amplitudes for intermediate states = 0
    int current_state = init_state;
    int next_state, l;
	int* colours = new int[(L-1)*NQ];
	colours = fs_bs_RG(circuit);
    
    for (l=0 ; l<L ; l++) {
        float lR=1.f, lI=0.f;
        
        next_state = (l< L-1 ? ndxs[l] : final_state);
        
        // processar o layer (lR, lI)
        // get gate layer l
        TCircuitLayer *layer = &circuit->layers[l];
        
        
        
        
        
        // get the amplitude of this layer for the given current and next states
        // We are reevaluating all layers even though a single ndxs[] changes: INNEFICIENT
        layer_w (layer, l, current_state, next_state, lR, lI);
        
        layer_aR[l] = lR;
        layer_aI[l] = lI;

        current_state = next_state;

    }

    
    while (ndxs[0] < N) {  // THIS IS THE EXPONENTIAL LOOP

        // .... process path : init_state-> ndxs[0] -> ..... -> ndxs[L-2] -> final_state
        
        float pathR=1.f, pathI=0.f;
                
        for (l=0 ; l<L ; l++) {  // multiply the layers' amplitudes
            
            complex_multiply (pathR, pathI, layer_aR[l], layer_aI[l], pathR, pathI);

        }
        
        sumR += pathR;
        sumI += pathI;

        
        // iterate to next path
        int ll;
        bool stop_updating = false;
        for (ll=L-2 ; ll>=0 ; ll--) {
            ndxs[ll]++;
            if (ndxs[ll]==N && ll>0)  {  // this loop resets
                ndxs[ll] = 0;
            }
            else
                stop_updating = true;
            
            // recompute layers associated with the changed intermediate
            // ndxs[ll] changed for sure
            
            // update layer_aR,I [ll+1]
            float lR=1.f, lI=0.f;
            
            current_state = ndxs[ll] ;
            next_state = (ll == L-2 ? final_state : ndxs[ll+1] );
            
            // processar o layer (lR, lI)
            // get gate layer l
            TCircuitLayer *layer = &circuit->layers[ll+1];
            
            // get the amplitude of this layer for the given current and next states
            // We are reevaluating all layers even though a single ndxs[] changes: INNEFICIENT
            layer_w (layer, ll+1, current_state, next_state, lR, lI);
            
            layer_aR[ll+1] = lR;
            layer_aI[ll+1] = lI;


            // update layer_aR,I [ll]
            lR=1.f;
            lI=0.f;
            
            current_state = (ll == 0 ? init_state : ndxs[ll-1] );
            next_state = ndxs[ll] ;
            
            // processar o layer (lR, lI)
            // get gate layer l
            layer = &circuit->layers[ll];
            
            // get the amplitude of this layer for the given current and next states
            // Wer are reevaluating all layers even though a single ndxs[] changes: INNEFICIENT
            layer_w (layer, ll, current_state, next_state, lR, lI);
            
            layer_aR[ll] = lR;
            layer_aI[ll] = lI;

            
            if (stop_updating) break;
        }

    }

    aR = sumR;
    aI = sumI;
    
}
