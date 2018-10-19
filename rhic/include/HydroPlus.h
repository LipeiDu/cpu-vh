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
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"

#define NUMBER_SLOW_MODES 3

typedef struct
{
    PRECISION *phiQ[NUMBER_SLOW_MODES];
} SLOW_MODES;

extern SLOW_MODES *PhiQ, *PhiQp, *PhiQS; // Slow modes: updated, previous, intermediate values
extern SLOW_MODES *eqPhiQ, *eqPhiQp, *eqPhiQS;  // Slow modes at equilibrium: updated, previous, intermediate values

extern PRECISION *Qvec; // Q vectors of slow modes

void allocateHostMemorySlowModes(int len);

PRECISION Cp(PRECISION e, PRECISION rhob);
PRECISION xi(PRECISION e, PRECISION rhob);//correlation length
PRECISION f2(PRECISION x); // Universal function f2, Eq.(89)

PRECISION equilibriumPhi0(PRECISION e, PRECISION rhob); // magnitude of the fluctuation at equilibirium, Eq.(90)
PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION Q); // cf. Eq.(89)

PRECISION relaxationCoefficientPhiQ(PRECISION e, PRECISION rhob, PRECISION Q);

void setInitialConditionSlowModes(void * latticeParams, void * hydroParams);

PRECISION entropySlowModes(const PRECISION * const __restrict__ equilibriumPhiQ, const PRECISION * const __restrict__ PhiQ, const PRECISION * const __restrict__ Q);

void setSlowModesSourceTerms(PRECISION * const __restrict__ PhiQRHS, const PRECISION * const __restrict__ equilibriumPhiQ, const PRECISION * const __restrict__ PhiQ, const PRECISION * const __restrict__ Qvec, PRECISION e, PRECISION rhob, PRECISION ut, PRECISION dkvk);

void loadSlowModesSourceTerms(const SLOW_MODES * const __restrict__ equilibriumSlowModes, const SLOW_MODES * const __restrict__ SlowModes, PRECISION * const __restrict__ S, const PRECISION * const __restrict__ Qvec, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ rhobvec, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz);


#endif /* HydroPlus_h */
