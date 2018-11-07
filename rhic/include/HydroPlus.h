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

PRECISION xi(PRECISION e, PRECISION rhob);

PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION s, PRECISION Q);

PRECISION relaxationCoefficientPhi(PRECISION rhob, PRECISION s, PRECISION T, PRECISION corrL2);

PRECISION relaxationCoefficientPhiQ(PRECISION gammaPhi, PRECISION corrL2, PRECISION Q);

void setInitialConditionSlowModes(void * latticeParams, void * hydroParams);

void getPrimaryVariablesFromSlowModes(PRECISION * const __restrict__ p, PRECISION * const __restrict__ T, PRECISION * const __restrict__ alphaB, const PRECISION * const __restrict__ equiPhiQ, const PRECISION * const __restrict__ PhiQ, PRECISION eIn, PRECISION rhobIn, PRECISION pIn, PRECISION TIn, PRECISION alphaBIn);

#endif /* HydroPlus_h */
