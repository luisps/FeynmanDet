//
//  gates.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 01/08/2023.
//

#ifndef gates_h
#define gates_h

#include "../States.h"

// math constants according to
// https://www.quantstart.com/articles/Mathematical-Constants-in-C/
#define _USE_MATH_DEFINES
#include <cmath>

const float M_SQRT2_OVER_2 = 0.707106781186548f;
const float MY_SQRT1_2 = 0.707106781186548f;

// Correspondence between gate names numbers and strings
// G1P0 - 1 qubit, no parameters
// 'id' - 0
// 'h'  - 1 B
// 'x'  - 2
// 'y'  - 3
// 'z'  - 4
// 's'  - 5 
// 't'  - 6
// G1P1 - 1 qubit, 1 parameter
// 'rx'  - 11 B
// 'ry'  - 12 B
// 'rz'  - 13
// 'p'  - 14
// G2P0 - 2 qubits, no parameters
// 'id2'  - 20    Used for errors
// 'cx'  - 21 B
// 'cz'  - 22
// G2P1 - 2 qubits, 1 parameter
// 'cp'  - 31

// G1P0 gates
inline void gate_id_w (int i_qb, int o_qb, float &wR, float &wI) {
    wR = (i_qb == o_qb ? 1.f : 0.f);
    wI = 0.f;
}
  
inline void gate_h_w (int i_qb, int o_qb, float &wR, float &wI) {
    wR = ((i_qb == 1)  && (o_qb == 1) ? -MY_SQRT1_2 : MY_SQRT1_2);
    wI = 0.f;
}

inline void gate_x_w (int i_qb, int o_qb, float &wR, float &wI) {
    wR = (i_qb == o_qb  ? 0.f : 1.f);
    wI = 0.f;
}

inline void gate_y_w (int i_qb, int o_qb, float &wR, float &wI) {
    wI = (i_qb == o_qb  ? 0.f : (i_qb == 0 ? 1.f : -1.f ));
    wR = 0.f;
}

inline void gate_z_w (int i_qb, int o_qb, float &wR, float &wI) {
    wR = (i_qb != o_qb  ? 0.f : (i_qb == 0 ? 1.f : -1.f ));
    wI = 0.f;
}

inline void gate_s_w (int i_qb, int o_qb, float &wR, float &wI) {
    wR = ((i_qb == 0) && (o_qb==0)  ? 1.f : 0.f );
    wI = ((i_qb == 1) && (o_qb==1)  ? 1.f : 0.f );
}

inline void gate_t_w (int i_qb, int o_qb, float &wR, float &wI) {
    wR = ((i_qb == o_qb)  ? (i_qb == 0 ? 1.f : M_SQRT2_OVER_2) : 0.f );
    wI = ((i_qb == o_qb)  ? (i_qb == 0 ? 0.f : M_SQRT2_OVER_2) : 0.f );
}

// G1P1 gates attached 2x2 unitary is used
inline void gate_g1p1_w (int i_qb, int o_qb,  float m[2][2][2], float &wR, float &wI) {
    //input selects the column
    // output selects the row
    wR = m[o_qb][i_qb][0];  // REAL
    wI = m[o_qb][i_qb][1];  // IMAG

}


// G2P0 gates
inline void gate_id2_w (int i_qbs, int o_qbs, float &wR, float &wI) {
    wR = (i_qbs == o_qbs ? 1.f : 0.f);
    wI = 0.f;
}

inline void gate_cx_w (int i_qbs, int o_qbs, float &wR, float &wI) {
    wR = (i_qbs < 2 && i_qbs == o_qbs ? 1.f : ((i_qbs==2 && o_qbs==3) || (i_qbs==3 && o_qbs ==2) ? 1.f : 0.f));
    wI = 0.f;
}

inline void gate_cz_w (int i_qbs, int o_qbs, float &wR, float &wI) {
    wR = (i_qbs < 3 && i_qbs == o_qbs ? 1.f : (i_qbs==3 && o_qbs==3 ? -1.f : 0.f));
    wI = 0.f;
}


// G2P1 gates attached 4x4 unitary is used
inline void gate_g2p1_w (int i_qbs, int o_qbs,  float m[4][4][2], float &wR, float &wI) {
    //input selects the column
    // output selects the row
    wR = m[o_qbs][i_qbs][0];  // REAL
    wI = m[o_qbs][i_qbs][1];  // IMAG

}


/* qubit handling functions for unsigned long long states */

inline int qb_value (int qb, StateT state) {
    const StateT uvalue = 1 <<qb;
    return ((state & uvalue) >> qb);
}

inline int qb_set (int qb, StateT state) {
    const StateT uvalue = (1ull <<qb);
    return (state | uvalue);
    
}

inline int qb_reset (int qb, StateT state) {
    const StateT uvalue = ~(1 <<qb);
    return state & uvalue;
    
}

#endif /* gates_h */
