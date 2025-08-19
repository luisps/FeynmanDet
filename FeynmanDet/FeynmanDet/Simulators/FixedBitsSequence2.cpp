//
//  FixedBitsSequence.cpp
//  IterateStates
//
//  Created by Luis Paulo Santos on 09/08/2025.
//

#include "FixedBitsSequence2.hpp"


// Assumes little endian
void FixedBitsSequence2::printBits(void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    int size = (nbits+7)/8;
    
    for (i = size-1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

// Generate all integers of 'n_bits' with bits fixed according to first value and non fixed according to nonfixed_bits
// This method returns the integer with index seq_ndx in that sequence
bool FixedBitsSequence2::generate_next_in_sequence(StateT &value) {

    if (seq_finished()) return false;
    
    StateT candidate = first_value;

    if (seq_ndx!=0) {
        for (int bit = 0; bit < nNonfixed_bits; bit++) {
            int bitval = (seq_ndx & (1U << bit));
            if (bitval !=0) { // '1' in seq_ndx, position bit
                candidate |= (1U << nonfixed_bits[bit]);
            }
        }
    }
    seq_ndx++;
    value = candidate;
    return true;
}




