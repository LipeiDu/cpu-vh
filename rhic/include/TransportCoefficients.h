//
//  TransportCoefficients.hpp
//  
//
//  Created by Lipei Du on 3/10/19.
//

#ifndef TransportCoefficients_h
#define TransportCoefficients_h

#include <stdio.h>

//Transport coefficients of the baryon evolution; Lipei
const PRECISION Cb = 4.0;

//Transport coefficients
const PRECISION delta_pipi = 1.33333;
const PRECISION tau_pipi = 0;//1.42857;
const PRECISION delta_PiPi = 0.666667;
const PRECISION lambda_piPi = 1.2;

PRECISION bulkViscosityToEntropyDensity(PRECISION T);


// baryon diffusion coefficients

void getBaryonDiffusionCoeffTable();

void baryonDiffusionCoeff(PRECISION T, PRECISION muB, PRECISION * const __restrict__ diffusionCoeff);

PRECISION baryonDiffusionCoefficient(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p);

PRECISION criticalBaryonDiffusionCoefficientAdscft(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq, PRECISION corrL);

PRECISION criticalBaryonDiffusionCoefficientPlus(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq);


#endif /* TransportCoefficients_hpp */
