//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#ifndef __SIMULATOR_RG__
#define __SIMULATOR_RG__

#include "../States.h"
#include "../BaseFunctionalty/circuit.h"

void simulate_RG_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI);

#endif
