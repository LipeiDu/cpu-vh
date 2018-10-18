//
//  HydroPlus.h
//  
//
//  Created by Lipei Du on 10/17/18.
//

// Equation indices from PRD 98 (2018) 036006

#ifndef HydroPlus_h
#define HydroPlus_h

#include "../include/DynamicalVariables.h"

#define NUMBER_SLOW_MODES 3

typedef struct
{
    PRECISION *qvec[NUMBER_SLOW_MODES];
} SLOW_MODES_Q;

typedef struct
{
    PRECISION *phiQ[NUMBER_SLOW_MODES];
} SLOW_MODES;

extern SLOW_MODES_Q *Qvec;
extern SLOW_MODES *phiQ;

//correlation length
PRECISION xi(PRECISION e, PRECISION rhob);

//fluctuations at equilibrium
PRECISION f2(PRECISION x); // Universal function f2, Eq.(89)
PRECISION equilibriumPhi0(PRECISION e, PRECISION rhob); // magnitude of the fluctuation at equilibirium, Eq.(90)
PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION Q); // cf. Eq.(89)

//relaxation of fluctuations
PRECISION relaxationCoefficientPhiQ(PRECISION e, PRECISION rhob, PRECISION Q);

//source terms of extra modes
void setSlowModesSourceTerms(PRECISION * const __restrict__ PhiQRHS, const PRECISION * const __restrict__ equilibriumPhiQ, const PRECISION * const __restrict__ PhiQ, const SLOW_MODES_Q * const __restrict__ Qvec, PRECISION e, PRECISION rhob, PRECISION ut, PRECISION dkvk);

void loadSlowModesSourceTerms(const SLOW_MODES * const __restrict__ equilibriumSlowModes, const SLOW_MODES * const __restrict__ SlowModes, PRECISION * const __restrict__ S, const SLOW_MODES_Q * const __restrict__ Qvec, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ rhobvec, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz);


#endif /* HydroPlus_h */
