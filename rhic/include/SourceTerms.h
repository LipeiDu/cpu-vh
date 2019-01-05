/*
 * SourceTerms.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef SOURCETERMS_H_
#define SOURCETERMS_H_

#include "../include/DynamicalVariables.h"
PRECISION baryonDiffusionCoefficient(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p);

PRECISION criticalBaryonDiffusionCoefficientAdscft(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq, PRECISION corrL);

PRECISION criticalBaryonDiffusionCoefficientPlus(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq);

void loadSourceTermsX(const PRECISION * const __restrict__ I, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dx);

void loadSourceTermsY(const PRECISION * const __restrict__ J, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dy);

void loadSourceTermsZ(const PRECISION * const __restrict__ K, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION t,
PRECISION d_dz);

void loadSourceTerms2(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp, PRECISION t, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ pvec, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_etabar, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz, const DYNAMICAL_SOURCE * const __restrict__ Source, const PRECISION * const __restrict__ rhobvec, const PRECISION * const __restrict__ alphaBvec, const PRECISION * const __restrict__ alphaBp, const PRECISION * const __restrict__ Tvec, PRECISION Tp, PRECISION seq, const SLOW_MODES *  const __restrict__ eqPhiQ);


#endif /* SOURCETERMS_H_ */
