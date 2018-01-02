/*
 * SourceTerms.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef SOURCETERMS_H_
#define SOURCETERMS_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

void loadSourceTerms(
const PRECISION * const __restrict__ I, const PRECISION * const __restrict__ J, const PRECISION * const __restrict__ K, 
const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S,
const PRECISION * const __restrict__ utvec, const PRECISION * const __restrict__ uxvec, 
const PRECISION * const __restrict__ uyvec, const PRECISION * const __restrict__ unvec,
PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp,
PRECISION t, PRECISION e, const PRECISION * const __restrict__ pvec,
int s
);
//=================================================================
void loadSourceTermsX(const PRECISION * const __restrict__ I, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dx
);

void loadSourceTermsY(const PRECISION * const __restrict__ J, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dy
);

void loadSourceTermsZ(const PRECISION * const __restrict__ K, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION t,
PRECISION d_dz
);

void loadSourceTerms2(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u,
PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp,
PRECISION t, const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ pvec,
int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_etabar, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz,
const PARTICLE_SOURCE * const __restrict__ Part, const PRECISION * const __restrict__ rhob//Part & rhob by Lipei
);


#endif /* SOURCETERMS_H_ */
