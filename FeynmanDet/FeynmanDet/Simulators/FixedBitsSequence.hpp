//
//  FixedBitsSequence.hpp
//  IterateStates
//
//  Created by Luis Paulo Santos on 09/08/2025.
//

#ifndef FixedBitsSequence_hpp
#define FixedBitsSequence_hpp

#include <stdio.h>

#include "../States.h"

class FixedBitsSequence {       // The class
private:
    StateT fixed_mask, first_value;
    int free_bits, n_bits;
    // Map from position in free_bits to actual bit position
    int free_positions[sizeof(StateT)];
    StateT total_in_sequence;

    void eval_fixed_mask (int const n_fixed_bits, int const fixed_bits[]) {
        fixed_mask = 0;
        for (int i=0; i < n_fixed_bits ; i++) {
            fixed_mask |= (1U << fixed_bits[i]);
        }
    }

    void eval_first_value (int const n_fixed_bits, int const fixed_bits[], int const fixed_values[]) {
        first_value = 0;
        for (int i=0; i < n_fixed_bits ; i++) {
            first_value |= (fixed_values[i] << fixed_bits[i]);
        }
    }

    StateT eval_total_in_sequence (void) {

        for (int i = 0; i < n_bits; ++i) {
            if ((fixed_mask & (1U << i)) == 0) {
                free_positions[free_bits++] = i;
            }
        }

        total_in_sequence = 1U << free_bits;
        return total_in_sequence;
    }
public:             // Access specifier
    FixedBitsSequence () {};
    StateT init_iterator (int const _n_bits,
                        int const n_fixed_bits,
                        int fixed_bits[],
                        int fixed_values[]) {
        
        n_bits = _n_bits;
        eval_fixed_mask(n_fixed_bits, fixed_bits);
        eval_first_value (n_fixed_bits, fixed_bits, fixed_values);
        return eval_total_in_sequence ();
    }
    void printBits(void const * const ptr);
    
    StateT generate_next_in_sequence(StateT & ndx);
    
    StateT get_first_value (void) {return first_value;};
};


#endif /* FixedBitsSequence_hpp */
