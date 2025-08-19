//
//  FixedBitsSequence.hpp
//  IterateStates
//
//  Created by Luis Paulo Santos on 09/08/2025.
//

#ifndef FixedBitsSequence2_hpp
#define FixedBitsSequence2_hpp

#include <stdio.h>

#include "../States.h"

class FixedBitsSequence2 {       // The class
private:
    StateT first_value;
    int nfixed_bits, nbits, nNonfixed_bits;
    StateT total_in_sequence;  // 2^nNonfixed_bits
    // which are the non fixed
    int nonfixed_bits[sizeof(StateT)];
    int seq_ndx;  // counter of the last generated value in the seq

    void eval_first_value (int const fixed_bits[], int const fixed_values[]) {
        first_value = 0;
        for (int i=0; i < nfixed_bits ; i++) {
            first_value |= (fixed_values[i] << fixed_bits[i]);
        }
    }

public:             // Access specifier
    FixedBitsSequence2 () {};
    StateT init_iterator (int const _nbits,
                        int const _nfixed_bits,
                        int fixed_bits[],
                        int fixed_values[]) {
        
        nbits = _nbits;
        nfixed_bits = _nfixed_bits;
        nNonfixed_bits = nbits - nfixed_bits;
        total_in_sequence = (1U << nNonfixed_bits);
        // evaluate firts value in allowed sequence
        // all non fixed bits are zero
        // the fixed bits will have their fixed value (daahhh)
        eval_first_value (fixed_bits, fixed_values);
        // compute which are the non fixed bits
        for (int i=0, nf=0 ; i < nbits ; i++) {
            bool found = false;
            for (int fix=0 ; !found && fix < nfixed_bits ; fix++) {
                if (i==fixed_bits[fix]) found = true;
            }
            if (!found) {  // it is non fixed
                nonfixed_bits[nf] = i;
                nf++;
            }
        }
        seq_ndx = 0;
        return first_value;
    }
    void printBits(void const * const ptr);
    
    // true if a valid number is generated
    // false if end of sequence reached
    bool generate_next_in_sequence(StateT &value);
    
    // to verify loop terminations
    bool seq_finished (void) {
        return (seq_ndx == total_in_sequence);
    }
    
    void reset(void) {seq_ndx=0;}
    StateT get_first_value (void) {return first_value;};
    StateT get_seq_length (void) {return total_in_sequence;};
};


#endif /* FixedBitsSequence2_hpp */
