//  Path_simulator
//
//  Created by Luis Paulo Santos and Vitoria Sousa on 04/12/2023.

#ifndef __SIMULATOR_PB__
#define __SIMULATOR_PB__

#include "../States.h"
#include "../BaseFunctionalty/circuit.h"

void simulate_PB_paths (TCircuit *circuit, StateT init_state, StateT final_state, float& aR, float& aI);

#endif
