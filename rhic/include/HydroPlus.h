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

//correlation length
PRECISION xi(PRECISION e, PRECISION rhob);

//fluctuations at equilibrium
PRECISION f2(PRECISION x); // Universal function f2, Eq.(89)
PRECISION equilibriumPhi0(PRECISION e, PRECISION rhob); // magnitude of the fluctuation at equilibirium, Eq.(90)
PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION Q); // cf. Eq.(89)

//relaxation of fluctuations
PRECISION relaxationCoefficientPhiQ(PRECISION e, PRECISION rhob, PRECISION Q);

//source terms of extra modes
void setPhiQsourceTerms(PRECISION * const __restrict__ PhiQRHS, PRECISION * const __restrict__ equilibriumPhiQRHS, const Slow_Modes * const __restrict__ phiQ,
                        PRECISION T, PRECISION t, PRECISION e, PRECISION p, PRECISION ut, PRECISION dkvk);


#endif /* HydroPlus_h */
