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

    if( generate_this_in_sequence(value, seq_ndx)) {
        seq_ndx++;
        return true;
    }
    else {
        return false;
    }
}

// Generate all integers of 'n_bits' with bits fixed according to first value and non fixed according to nonfixed_bits
// This method returns the integer with index this_seq_ndx in that sequence
bool FixedBitsSequence2::generate_this_in_sequence(StateT &value, StateT const this_seq_ndx) {

    if (seq_finished(this_seq_ndx)) return false;
    
    StateT candidate = first_value;

    if (this_seq_ndx!=0) {
        for (int bit = 0; bit < nNonfixed_bits; bit++) {
            int bitval = (this_seq_ndx & (1U << bit));
            if (bitval !=0) { // '1' in this_seq_ndx, position bit
                candidate |= (1U << nonfixed_bits[bit]);
            }
        }
    }
    value = candidate;
    return true;
}



