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

/*static int gate_output_1(int qbv, int name){
    int out=0;
    
    switch (name){
        case 0:
        case 4:
        case 5:
        case 6:
        case 13:
        { //id
            if(qbv==1){
                out=1;
            }
            break;
        }
            
        case 2:
        case 3:
        { //X //Y
            if(qbv==0){
                out=1;
            }
            break;
        }
        default:
            fprintf (stderr, "gate_output error: unknown or branching gate found!\n");
            break;
    } // switch
    return out;
}


// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
// gates 2 qubits
// G2P0 - 2 qubits, no parameters
// 'id2'  - 20    Used for errors
// 'cx'  - 21 B
// 'cz'  - 22
// G2P1 - 2 qubits, 1 parameter
// 'cp'  - 31

static void gate_output_2(int qbv_c, int qbv_t, int name, int* out){
    
    out[0]=0;
    out[1]=0;
    switch (name){
        case 20:
        case 22:
        case 31: {
            out[0]=qbv_c;
            out[1]=qbv_t;
            break;
        }
            
        case 21: { //CX
            if(qbv_c==1){
                out[0]=1;
                if(qbv_t==0){
                    out[1]=1;
                }
            } else {
                if(qbv_t==1){
                    out[1]=1;
                }
            }
            
            break;
        }
    }
}


static int qbv_given_previous_gate (TCircuit *circuit, int l, int n_qb, StateT previous_state) {
    int qbv;
        
    TCircuitLayer* layer = &circuit->layers[l];
    // Iterate over the gate types for each layer
    for (int GT = 0; GT < 4; GT++) {
        void* gates_ptr = layer->gates[GT];
        int num_gates = layer->num_type_gates[GT];
        switch (GT) {
            case 0: {
                TGate1P0* g1p0_ptr = (TGate1P0*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p0_ptr++) {
                    int qubit = g1p0_ptr->qubit; //qubit que participa na gate
                    if (qubit==n_qb) {
                        int name = g1p0_ptr->name;
                        int prev_qbv = qb_value(n_qb, previous_state);
                        qbv = gate_output_1(prev_qbv, name);
                        return qbv;
                    }
                }
                break;
            }
            case 1: {
                TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                    int qubit = g1p1_ptr->fdata.qubit; //qubit que participa na gate
                    if (qubit==n_qb) {
                        int name = g1p1_ptr->fdata.name;
                        int prev_qbv = qb_value(n_qb, previous_state);
                        qbv = gate_output_1(prev_qbv, name);
                        //printf ("gate name: %d ; input=%d; out=%d;\n",name, prev_qbv,qbv);
                       return qbv;
                    }
                }
                break;
            }
            case 2: {
                TGate2P0* g2p0_ptr = (TGate2P0*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g2p0_ptr++) {
                    int c_qubit = g2p0_ptr->c_qubit;
                    int t_qubit = g2p0_ptr->t_qubit;
                    if (c_qubit==n_qb) {
                        //printf ("previous_state: %llu ; c_ndx: %d ; t_ndx: %d\n",previous_state, c_qubit, t_qubit);
                        int name = g2p0_ptr->name;
                        int prev_qbv_c = qb_value(c_qubit, previous_state);
                        int prev_qbv_t = qb_value(t_qubit, previous_state);
                        int out[2];
                        gate_output_2(prev_qbv_c, prev_qbv_t, name, out);
                        qbv = out[0];
                        //printf ("gate name: %d ; input: c=%d, t=%d; out: c=%d, t=%d;\n",name, prev_qbv_c,prev_qbv_t,out[0],out[1]);
                        return qbv;
                    }
                    if (t_qubit==n_qb) {
                        int name = g2p0_ptr->name;
                        int prev_qbv_c = qb_value(c_qubit, previous_state);
                        int prev_qbv_t = qb_value(t_qubit, previous_state);
                        int out[2];
                        gate_output_2(prev_qbv_c, prev_qbv_t, name, out);
                        qbv = out[1];
                        return qbv;
                    }
                }
                break;
            }
            case 3: {
                TGate2P1* g2p1_ptr = (TGate2P1*)gates_ptr;
                for (int g = 0; g < num_gates; g++, g2p1_ptr++) {
                    int c_qubit = g2p1_ptr->fdata.c_qubit; //qubit que participa na gate
                    int t_qubit = g2p1_ptr->fdata.t_qubit;
                    if (c_qubit==n_qb) {
                        int name = g2p1_ptr->fdata.name;
                        int prev_qbv_c = qb_value(c_qubit, previous_state);
                        int prev_qbv_t = qb_value(t_qubit, previous_state);
                        int out[2];
                        gate_output_2(prev_qbv_c, prev_qbv_t, name, out);
                        qbv = out[0];
                        return qbv;
                    }
                    if (t_qubit==n_qb) {
                        int name = g2p1_ptr->fdata.name;
                        int prev_qbv_c = qb_value(c_qubit, previous_state);
                        int prev_qbv_t = qb_value(t_qubit, previous_state);
                        int out[2];
                        gate_output_2(prev_qbv_c, prev_qbv_t, name, out);
                        qbv = out[1];
                        return qbv;
                    }
                }
                break;
            }
        } // switch
    } // for GT
    return 0;
}*/

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

void simulate_PB_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI) {

    const int L = circuit->size->num_layers;
    StateT *ndxs = new StateT[L-1];
    const int NQ=circuit->size->num_qubits;
    const StateT N = 1 << NQ;
    StateT path_counter=0, path_NZ_counter=0;

    float* wR=new float[L-1];
    float* wI=new float[L-1];

    /*for (StateT ant=0 ; ant <7 ; ant++) {
        int blue_qbv = qbv_given_previous_gate (circuit, 1, 0, ant);
        fprintf (stderr, "input=%d ; output=%d\n", qb_value(0, ant), blue_qbv);
    }
    return;*/

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

    fs_PB(circuit, colours);
    
    print_PB (colours, NQ, L);
    int start_layer=0;

    // Simulation starts

    // all intermediate layers indexes to 0
    for (int i=0 ; i<L-1 ; i++) ndxs[i]=0 ;

    // main simulation loop
    while (ndxs[0] < N) {

        
        //fprintf (stderr, "ndxs= ");
        //for (int lll=0; lll<L-1 ; lll++)
        //    fprintf (stderr, "\t%05llu", ndxs[lll]);
        //fprintf (stderr, "\n");

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

        } // for all layers
        
        if  (!zero_weight_layer) {
            sumR += pathR;
            sumI += pathI;
            path_NZ_counter++;
            // DEBUG
            /*
            printf ("Non zero path: ");
            for (int lll=0 ; lll<L-1 ; lll++) {
                printf ("%llu ", ndxs[lll]);
            }
            printf ("= %e + i %e\n", pathR, pathI);*/
        }
        path_counter++;

        // compute next path
        // updating ndxs[]
    
        // compute next path skipping invalid GREENs and fixed BLUEs
        int ll;
        bool invalid_state_green = true;
        for (ll=((zero_weight_layer && l<(L-1))? l : L-2); ll>=0 && invalid_state_green ; ll--) {

            //printf ("change ndxs[%d]=%llu\n", ll, ndxs[ll]);
            // get the state from the previous layer. Will need it
            StateT const prev_state = (ll==0 ? init_state : ndxs[ll-1]);

            invalid_state_green = true;
            while(invalid_state_green) {
                ndxs[ll]++;
                if (ll==0)  fprintf (stderr, "ndxs[0] = %llu \n", ndxs[0]);

                if (ndxs[ll]==N && ll>0)  { // this layer overflows
                    ndxs[ll] = 0;
                    break;        // break only from inner loop
                }
                else if (ndxs[0]==N && ll==0)  { // simulation finished
                    invalid_state_green = false; // terminate outer loop
                    break;   // terminate inner loop
                }
                start_layer=ll;
                //printf ("changed ndxs[%d]=%llu (ndxs[0] = %llu) \n", ll, ndxs[ll], ndxs[0]);
                invalid_state_green = false;
                // verify whether this ndxs complies with the colouring
                for (int i=0; i<NQ && !invalid_state_green ; i++){
                    int const next_state_CL = colours[i+ll*NQ];
                    
                    if (next_state_CL== GREEN0 || next_state_CL== GREEN1) {
                        if (next_state_CL != qb_value(i,ndxs[ll])) {
                            invalid_state_green=true; // break from inner loop
                            ndxs[ll] |= ((1 << i)- 1) ; // skip all  intermediate non valid states
                        }
                    }
                    else if (next_state_CL == BLUE_I || next_state_CL == BLUE_X) {
                        int const pred_qb_value = (next_state_CL == BLUE_I ? qb_value(i,prev_state) : !qb_value(i,prev_state));
                        if (pred_qb_value != qb_value(i,ndxs[ll])) {
                            invalid_state_green=true;
                            ndxs[ll] |= ((1 << i)- 1) ; // skip all  intermediate non valid states
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
                        if (pred_qb_value != qb_value(i,ndxs[ll])) {
                            invalid_state_green=true;
                            ndxs[ll] |= ((1 << i)- 1) ; // skip all  intermediate non valid states
                        }
                    }

                } // for to verify qubits and green
            }  // while (invalid_state_green)
            //printf ("END FOR LOOP ndxs[%d]=%llu\n", ll, ndxs[ll]);
        }    // for backward change layers ndxs
        
    } // main simulation loop (while)
    aR = sumR;
    aI = sumI;

    printf ("< %llu | U | %llu > = %.6f + i %.6f, p=%.6f\n", final_state, init_state, aR, aI, aR*aR+aI*aI);
    printf ("%llu paths, %llu non zero\n", path_counter, path_NZ_counter);
}
