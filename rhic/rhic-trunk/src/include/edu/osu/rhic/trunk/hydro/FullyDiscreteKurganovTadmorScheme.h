/*
 * FullyDiscreteKurganovTadmorScheme.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef FULLYDISCRETEKURGANOVTADMORSCHEME_H_
#define FULLYDISCRETEKURGANOVTADMORSCHEME_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

void rungeKutta2(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q, 
void * latticeParams, void * hydroParams
);
void eulerStepKernelSourcePart(PRECISION t, PRECISION dt, const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,int ncx, int ncy, int ncz);// by Lipei

#endif /* FULLYDISCRETEKURGANOVTADMORSCHEME_H_ */
