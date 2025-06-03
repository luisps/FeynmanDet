//
//  simulate_all.hpp
//  Path_simulator
//
//  Created by Luis Paulo Santos on 03/11/2023.
//

#ifndef simulate_all_hpp
#define simulate_all_hpp

#include <stdio.h>
#include "circuit.h"

void simulate_all_paths (TCircuit *circuit,
                         int init_state, int final_state,
                         float& aR, float& aI);

void simulate_all_paths_opt (TCircuit *circuit,
                         int init_state, int final_state,
                         float& aR, float& aI);

#endif /* simulate_all_hpp */
