//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#ifndef __SIMULATOR_RG_NEW__
#define __SIMULATOR_RG_NEW__

#include "../States.h"
#include "../BaseFunctionalty/circuit.h"
#include "FixedBitsSequence2.hpp"

void simulate_RG_paths_new_errors (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI);

#endif
