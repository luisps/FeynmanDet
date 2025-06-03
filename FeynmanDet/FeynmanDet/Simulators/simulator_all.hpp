//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#ifndef __SIMULATOR_ALL__
#define __SIMULATOR_ALL__

#include "../States.h"
#include "../BaseFunctionalty/circuit.h"

void simulate_all_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI);

#endif
