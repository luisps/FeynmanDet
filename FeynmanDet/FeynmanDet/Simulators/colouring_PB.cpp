//  Path_simulator
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.
//
#include "colouring_PB.hpp"

#include <math.h>
#include <string.h>
#include <stdio.h>

// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
//gates 1 qubit

// Correspondence between gate names numbers and strings
// G1P0 - 1 qubit, no parameters
// 'id' - 0 BLUE_I
// 'h'  - 1 Branching
// 'x'  - 2 BLUE_X
// 'y'  - 3 BLUE_X
// 'z'  - 4 BLUE_I
// 's'  - 5 BLUE_I
// 't'  - 6 BLUE_I
// G1P1 - 1 qubit, 1 parameter
// 'rx'  - 11 Branching
// 'ry'  - 12 Branching
// 'rz'  - 13 BLUE_I
// 'p'  - 14  BLUE_I


// pega no input estado inicial e gate -> d� o estado final (para ver se � 1 ou 0)
// gates 2 qubits
// G2P0 - 2 qubits, no parameters
// 'id2'  - 20    Used for errors
// 'cx'  - 21 control: BLUE_I; target BLUE_CX
// 'cz'  - 22 control: BLUE_I; target BLUE_I
// G2P1 - 2 qubits, 1 parameter
// 'cp'  - 31 BLUE_I



void fs_PB(TCircuit *circuit, colourT* states_colour){
	
//Layers e qubits    
    const int L = circuit->size->num_layers;
    const int NQ = circuit->size->num_qubits;
    
    for (int l=0 ; l<L ; l++) {
        TCircuitLayer* layer = &circuit->layers[l];
        // Iterate over the gate types for each layer
        for (int GT = 0; GT < 4; GT++) {
            void* gates_ptr = layer->gates[GT];
            int num_gates = layer->num_type_gates[GT];
            switch (GT) {
                    
                case 0: {
                    TGate1P0* g1p0_ptr = (TGate1P0*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g1p0_ptr++) {
                        int name = g1p0_ptr->name;
                        int qubit = g1p0_ptr->qubit; //qubit que participa na gate
                        colourT const colour= states_colour[l*NQ+qubit];
                        if (colour!=RED) continue;
                        // we are now dealing with a RED that can become
                        // PINK or BLUE
                        if (name==1) {  // Hadamard - branching
                            states_colour[l*NQ+qubit] = PINK;
                        } else if (name==2 || name==3) { // X or Y
                            //states_colour[l*NQ+qubit] = BLUE;
                            states_colour[l*NQ+qubit] = BLUE_X;
                        } else {
                            //states_colour[l*NQ+qubit] = BLUE;
                            states_colour[l*NQ+qubit] = BLUE_I;
                        }
                    }
                    break;
                }
                case 1: {
                    TGate1P1* g1p1_ptr = (TGate1P1*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g1p1_ptr++) {
                        int name = g1p1_ptr->fdata.name;
                        int qubit = g1p1_ptr->fdata.qubit; //qubit que participa na gate
                        colourT const colour= states_colour[l*NQ+qubit];
                        if (colour!=RED) continue;
                        // we are now dealing with a RED that can become
                        // PINK or BLUE
                        if (name==11 || name==12) {  // RX,RY - branching
                            states_colour[l*NQ+qubit] = PINK;
                        } else {
                            //states_colour[l*NQ+qubit] = BLUE;
                            states_colour[l*NQ+qubit] = BLUE_I;
                        }
                    }
                    break;
                }
                case 2: {
                    TGate2P0* g2p0_ptr = (TGate2P0*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g2p0_ptr++) {
                        int name = g2p0_ptr->name;
                        int c_qubit = g2p0_ptr->c_qubit;
                        int t_qubit = g2p0_ptr->t_qubit;
                        
                        colourT const c_colour= states_colour[l*NQ+c_qubit];
                        colourT const t_colour= states_colour[l*NQ+t_qubit];
                        // PINK or BLUE. No PINKs in this case
                        // control
                        states_colour[l*NQ+c_qubit] = (c_colour==RED ? BLUE_I : c_colour);
                        // is this a CX or CZ ?
                        if (name==21) { // CX
                            if (t_colour != RED) { // this is green
                                states_colour[l*NQ+t_qubit] = t_colour;
                            } else {  // RED becomes BLUE_CX
                                // encode the index of the control qubit on the MS bytes
                                states_colour[l*NQ+t_qubit] = BLUE_CX | (c_qubit << 8);
                                
                            }
                        } else { // CZ
                            states_colour[l*NQ+t_qubit] = (t_colour==RED ? BLUE_I : t_colour);
                        }
                    }
                    break;
                }
                case 3: {
                    TGate2P1* g2p1_ptr = (TGate2P1*)gates_ptr;
                    for (int g = 0; g < num_gates; g++, g2p1_ptr++) {
                        // necessarily CP
                        //int name = g2p1_ptr->fdata.name;
                        int c_qubit = g2p1_ptr->fdata.c_qubit; //qubit que participa na gate
                        int t_qubit = g2p1_ptr->fdata.t_qubit;
                        colourT const c_colour= states_colour[l*NQ+c_qubit];
                        colourT const t_colour= states_colour[l*NQ+t_qubit];
                        // PINK or BLUE. No PINKs in this case
                        states_colour[l*NQ+c_qubit] = (c_colour==RED ? BLUE_I : c_colour);
                        states_colour[l*NQ+t_qubit] = (t_colour==RED ? BLUE_I : t_colour);
                    }
                    break;
                }
            } // switch
        } // for GT
    } // for  l (layers)
}
