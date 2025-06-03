//
//  layer.hpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 18/08/2023.
//

#ifndef layer_hpp
#define layer_hpp

#include <stdio.h>

#include "circuit.h"
#include "gates.h"
#include "complex.h"

void layer_w (TCircuitLayer *layer, int l,
              int current_state, int next_state, float &wR, float &wI);


#endif /* layer_hpp */
