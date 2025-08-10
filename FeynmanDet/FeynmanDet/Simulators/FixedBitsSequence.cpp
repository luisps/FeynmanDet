//
//  FixedBitsSequence.cpp
//  IterateStates
//
//  Created by Luis Paulo Santos on 09/08/2025.
//

#include "FixedBitsSequence.hpp"


// Assumes little endian
void FixedBitsSequence::printBits(void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    int size = (n_bits+7)/8;
    
    for (i = size-1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

// Generate all integers of 'n_bits' with bits fixed according to 'fixed_mask' and 'fixed_value'
StateT FixedBitsSequence::generate_next_in_sequence(StateT &ndx) {

    StateT candidate = first_value;

    if (ndx!=0) {
        for (int j = 0; j < free_bits; ++j) {
            if (ndx & (1U << j)) {
                candidate |= (1U << free_positions[j]);
            } else {
                candidate &= ~(1U << free_positions[j]);
            }
        }
    }
    ndx++;
    return candidate;
}




