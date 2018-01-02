/*
 * EnergyMomentumTensor.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef ENERGYMOMENTUMTENSOR_H_
#define ENERGYMOMENTUMTENSOR_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

void getInferredVariables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION ePrev,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, 
PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un,
PRECISION rhobPrev, PRECISION * const __restrict__ rhob//rhob by Lipei
);

void setInferredVariablesKernel(const CONSERVED_VARIABLES * const __restrict__ q, 
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, 
PRECISION t, void * latticeParams, PRECISION * const __restrict__ rhob//rhob by Lipei
);

PRECISION Ttt(PRECISION e, PRECISION p, PRECISION ut, PRECISION pitt);
PRECISION Ttx(PRECISION e, PRECISION p, PRECISION ut, PRECISION ux, PRECISION pitx);
PRECISION Tty(PRECISION e, PRECISION p, PRECISION ut, PRECISION uy, PRECISION pity);
PRECISION Ttn(PRECISION e, PRECISION p, PRECISION ut, PRECISION un, PRECISION pitn);
PRECISION Txx(PRECISION e, PRECISION p, PRECISION ux, PRECISION pixx);
PRECISION Txy(PRECISION e, PRECISION p, PRECISION ux, PRECISION uy, PRECISION pixy);
PRECISION Txn(PRECISION e, PRECISION p, PRECISION ux, PRECISION un, PRECISION pixn);
PRECISION Tyy(PRECISION e, PRECISION p, PRECISION uy, PRECISION piyy);
PRECISION Tyn(PRECISION e, PRECISION p, PRECISION uy, PRECISION un, PRECISION piyn);
PRECISION Tnn(PRECISION e, PRECISION p, PRECISION un, PRECISION pinn, PRECISION t);

PRECISION Nbt(PRECISION rhob, PRECISION ut, PRECISION nbt);//Lipei
PRECISION Nbx(PRECISION rhob, PRECISION ux, PRECISION nbx);//Lipei
PRECISION Nby(PRECISION rhob, PRECISION uy, PRECISION nby);//Lipei
PRECISION Nbn(PRECISION rhob, PRECISION un, PRECISION nbn);//Lipei

#endif /* ENERGYMOMENTUMTENSOR_H_ */
