//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#ifndef __SIMULATOR_PB__
#define __SIMULATOR_PB__

#include "../States.h"
#include "../BaseFunctionalty/circuit.h"

#include "colouring_PB.hpp"

void simulate_PB_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI, char *NZ_paths_filename=NULL);

int validate_PB (StateT const state, int const *const colours, int const NQ);

#endif
