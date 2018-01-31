/*
 * FullyDiscreteKurganovTadmorScheme.cpp
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <stdlib.h>
#include <stdio.h> // for printf
#include <math.h>
#include <iostream>// by Lipei
#include <fstream>// by Lipei
#include <iomanip>//by Lipei


#include "edu/osu/rhic/trunk/hydro/FullyDiscreteKurganovTadmorScheme.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/core/muscl/SemiDiscreteKurganovTadmorScheme.h"
#include "edu/osu/rhic/core/muscl/HalfSiteExtrapolation.h"
#include "edu/osu//rhic/trunk/hydro/FluxFunctions.h"
#include "edu/osu//rhic/trunk/hydro/SpectralRadius.h"
#include "edu/osu/rhic/trunk/hydro/SourceTerms.h"
#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"

#include "edu/osu/rhic/core/util/FiniteDifference.h" //temp

using namespace std;//Lipei

/**************************************************************************************************************************************************/
void setNeighborCellsJK2(const PRECISION * const __restrict__ in, PRECISION * const __restrict__ out,
int s, int ptr, int smm, int sm, int sp, int spp
) {
	PRECISION data_ns= in[s];
	*(out + ptr    ) = in[smm];
	*(out + ptr + 1) = in[sm];
	*(out + ptr + 2) = data_ns;
	*(out + ptr + 3) = in[sp];
	*(out + ptr + 4) = in[spp];
}

void eulerStepKernelSource(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ ep, const PRECISION * const __restrict__ p,
const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dx, PRECISION dy, PRECISION dz, PRECISION etabar,
const PRECISION * const __restrict__ rhob, const PRECISION * const __restrict__ rhobp//ep, rhob, rhobp by Lipei
) {
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                
				PRECISION Q[ALL_NUMBER_CONSERVED_VARIABLES];
				PRECISION S[ALL_NUMBER_CONSERVED_VARIABLES];

				Q[0] = currrentVars->ttt[s];
				Q[1] = currrentVars->ttx[s];
				Q[2] = currrentVars->tty[s];
				Q[3] = currrentVars->ttn[s];
#ifdef PIMUNU
				Q[4] = currrentVars->pitt[s];
				Q[5] = currrentVars->pitx[s];
				Q[6] = currrentVars->pity[s];
				Q[7] = currrentVars->pitn[s];
				Q[8] = currrentVars->pixx[s];
				Q[9] = currrentVars->pixy[s];
				Q[10] = currrentVars->pixn[s];
				Q[11] = currrentVars->piyy[s];
				Q[12] = currrentVars->piyn[s];
				Q[13] = currrentVars->pinn[s];
#endif
#ifdef PI
				Q[14] = currrentVars->Pi[s];
#endif
#ifdef NBMU
                Q[NUMBER_CONSERVED_VARIABLES] = currrentVars->Nbt[s];
#endif
#ifdef VMU
                Q[NUMBER_CONSERVED_VARIABLES+1] = currrentVars->nbt[s];
                Q[NUMBER_CONSERVED_VARIABLES+2] = currrentVars->nbx[s];
                Q[NUMBER_CONSERVED_VARIABLES+3] = currrentVars->nby[s];
                Q[NUMBER_CONSERVED_VARIABLES+4] = currrentVars->nbn[s];
#endif

				loadSourceTerms2(Q, S, u, up->ut[s], up->ux[s], up->uy[s], up->un[s], t, e, ep[s], p, s, ncx, ncy, ncz, etabar, dt, dx, dy, dz, Source, rhob, rhobp[s]);
                
				PRECISION result[ALL_NUMBER_CONSERVED_VARIABLES];
				for (unsigned int n = 0; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = *(Q+n) + dt * ( *(S+n) );
				}

				updatedVars->ttt[s] = result[0];
				updatedVars->ttx[s] = result[1];
				updatedVars->tty[s] = result[2];
				updatedVars->ttn[s] = result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] = result[4];
				updatedVars->pitx[s] = result[5];
				updatedVars->pity[s] = result[6];
				updatedVars->pitn[s] = result[7];
				updatedVars->pixx[s] = result[8];
				updatedVars->pixy[s] = result[9];
				updatedVars->pixn[s] = result[10];
				updatedVars->piyy[s] = result[11];
				updatedVars->piyn[s] = result[12];
				updatedVars->pinn[s] = result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] = result[14];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] = result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] = result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] = result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] = result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] = result[NUMBER_CONSERVED_VARIABLES+4];
#endif

			}
		}
	}
}


void eulerStepKernelX(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dx, const PRECISION * const __restrict__ rhob
) {
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				PRECISION I[5 * ALL_NUMBER_CONSERVED_VARIABLES];
				PRECISION H[ALL_NUMBER_CONSERVED_VARIABLES];

				// calculate neighbor cell indices;
				int sim  = s-1;
				int simm = sim-1;
				int sip  = s+1;
				int sipp = sip+1;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,I,s,ptr,simm,sim,sip,sipp);  ptr+=5;
#endif
#ifdef NBMU
                setNeighborCellsJK2(currrentVars->Nbt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#endif
#ifdef VMU
                setNeighborCellsJK2(currrentVars->nbt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nby,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbn,I,s,ptr,simm,sim,sip,sipp);
#endif
                
				PRECISION result[ALL_NUMBER_CONSERVED_VARIABLES];
				flux(I, H, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusX, &Fx, t, e[s], rhob[s]);
				for (unsigned int n = 0; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = - *(H+n);
				}
				flux(I, H, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusX, &Fx, t, e[s], rhob[s]);
				for (unsigned int n = 0; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) += *(H+n);
					*(result+n) /= dx;
				}

                //===================================================
                // Conserved currents
                //===================================================
#ifndef IDEAL
                //------------Non-ideal case------------
                
				loadSourceTermsX(I, H, u, s, dx);
                
                //for energy-momentum tensor
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) += *(H+n);
					*(result+n) *= dt;
				}

                //for baryon current
#ifdef NBMU
                *(result+NUMBER_CONSERVED_VARIABLES) += *(H+NUMBER_CONSERVED_VARIABLES);
                *(result+NUMBER_CONSERVED_VARIABLES) *= dt;
#endif

#else
                //------------Ideal case------------
                
                //for energy-momentum tensor
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) *= dt;
				}
#ifdef NBMU
                //for baryon current
                *(result+NUMBER_CONSERVED_VARIABLES) *= dt;
#endif
#endif

                //===================================================
                // Dissipative currents
                //===================================================
                
#ifdef PIMUNU
				for (unsigned int n = 4; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) *= dt;
                }
#endif
#ifdef VMU
                for (unsigned int n = NUMBER_CONSERVED_VARIABLES+1; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
                    *(result+n) *= dt;
                }
#endif
                
                
				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[14];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] += result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] += result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] += result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] += result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] += result[NUMBER_CONSERVED_VARIABLES+4];
#endif
			}
		}
	}
}

void eulerStepKernelY(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dy, const PRECISION * const __restrict__ rhob
) {
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				PRECISION J[5* ALL_NUMBER_CONSERVED_VARIABLES];
				PRECISION H[ALL_NUMBER_CONSERVED_VARIABLES];

				// calculate neighbor cell indices;
				int sjm = s-ncx;
				int sjmm = sjm-ncx;
				int sjp = s+ncx;
				int sjpp = sjp+ncx;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,J,s,ptr,sjmm,sjm,sjp,sjpp);  ptr+=5;
#endif
#ifdef NBMU
                setNeighborCellsJK2(currrentVars->Nbt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#endif
#ifdef VMU
                setNeighborCellsJK2(currrentVars->nbt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nby,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbn,J,s,ptr,sjmm,sjm,sjp,sjpp);
#endif
                
                
				PRECISION result[ALL_NUMBER_CONSERVED_VARIABLES];
				flux(J, H, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusY, &Fy, t, e[s], rhob[s]);
				for (unsigned int n = 0; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = - *(H+n);
				}
				flux(J, H, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusY, &Fy, t, e[s], rhob[s]);
				for (unsigned int n = 0; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) += *(H+n);
					*(result+n) /= dy;
				}
#ifndef IDEAL
				loadSourceTermsY(J, H, u, s, dy);
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) += *(H+n);
					*(result+n) *= dt;
				}
#ifdef NBMU
                *(result+NUMBER_CONSERVED_VARIABLES) += *(H+NUMBER_CONSERVED_VARIABLES);
                *(result+NUMBER_CONSERVED_VARIABLES) *= dt;
#endif
#else
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) *= dt;
				}
#ifdef NBMU
                *(result+NUMBER_CONSERVED_VARIABLES) *= dt;
#endif
#endif

                
#ifdef PIMUNU
				for (unsigned int n = 4; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) *= dt;
				}
#endif
#ifdef VMU
                for (unsigned int n = NUMBER_CONSERVED_VARIABLES+1; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
                    *(result+n) *= dt;
                }
#endif

                
                
				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[14];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] += result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] += result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] += result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] += result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] += result[NUMBER_CONSERVED_VARIABLES+4];
#endif
			}
		}
	}
}

void eulerStepKernelZ(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dz, const PRECISION * const __restrict__ rhob
) {
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				PRECISION K[5 * ALL_NUMBER_CONSERVED_VARIABLES];
				PRECISION H[ALL_NUMBER_CONSERVED_VARIABLES];

				// calculate neighbor cell indices;
				int stride = ncx * ncy;
				int skm = s-stride;
				int skmm = skm-stride;
				int skp = s+stride;
				int skpp = skp+stride;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef NBMU
                setNeighborCellsJK2(currrentVars->Nbt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef VMU
                setNeighborCellsJK2(currrentVars->nbt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nby,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbn,K,s,ptr,skmm,skm,skp,skpp);
#endif
                
				PRECISION result[ALL_NUMBER_CONSERVED_VARIABLES];
				flux(K, H, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusZ, &Fz, t, e[s], rhob[s]);
				for (unsigned int n = 0; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = -*(H+n);
				}
				flux(K, H, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusZ, &Fz, t, e[s], rhob[s]);
				for (unsigned int n = 0; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) += *(H+n);
					*(result+n) /= dz;
				}
                
#ifndef IDEAL
				loadSourceTermsZ(K, H, u, s, t, dz);
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) += *(H+n);
					*(result+n) *= dt;
				}
#ifdef NBMU
                *(result+NUMBER_CONSERVED_VARIABLES) += *(H+NUMBER_CONSERVED_VARIABLES);
                *(result+NUMBER_CONSERVED_VARIABLES) *= dt;
#endif

#else
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) *= dt;
				}

#ifdef NBMU
                *(result+NUMBER_CONSERVED_VARIABLES) *= dt;
#endif
#endif
        
                
#ifdef PIMUNU
				for (unsigned int n = 4; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) *= dt;
				}
#endif
#ifdef VMU
                for (unsigned int n = NUMBER_CONSERVED_VARIABLES+1; n < ALL_NUMBER_CONSERVED_VARIABLES; ++n) {
                    *(result+n) *= dt;
                }
#endif

                
                
				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[14];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] += result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] += result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] += result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] += result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] += result[NUMBER_CONSERVED_VARIABLES+4];
#endif
			}
		}
	}
}

/**************************************************************************************************************************************************/

void convexCombinationEulerStepKernel(const CONSERVED_VARIABLES * const __restrict__ q, CONSERVED_VARIABLES * const __restrict__ Q,
int ncx, int ncy, int ncz
) {
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				Q->ttt[s] += q->ttt[s];
				Q->ttt[s] /= 2;
				Q->ttx[s] += q->ttx[s];
				Q->ttx[s] /= 2;
				Q->tty[s] += q->tty[s];
				Q->tty[s] /= 2;
				Q->ttn[s] += q->ttn[s];
				Q->ttn[s] /= 2;
#ifdef PIMUNU
				Q->pitt[s] += q->pitt[s];
				Q->pitt[s] /= 2;
				Q->pitx[s] += q->pitx[s];
				Q->pitx[s] /= 2;
				Q->pity[s] += q->pity[s];
				Q->pity[s] /= 2;
				Q->pitn[s] += q->pitn[s];
				Q->pitn[s] /= 2;
				Q->pixx[s] += q->pixx[s];
				Q->pixx[s] /= 2;
				Q->pixy[s] += q->pixy[s];
				Q->pixy[s] /= 2;
				Q->pixn[s] += q->pixn[s];
				Q->pixn[s] /= 2;
				Q->piyy[s] += q->piyy[s];
				Q->piyy[s] /= 2;
				Q->piyn[s] += q->piyn[s];
				Q->piyn[s] /= 2;
				Q->pinn[s] += q->pinn[s];
				Q->pinn[s] /= 2;
#endif
#ifdef PI
                Q->Pi[s]   += q->Pi[s];
                Q->Pi[s]   /= 2;
#endif
#ifdef NBMU
                Q->Nbt[s]  += q->Nbt[s];
                Q->Nbt[s]  /= 2;
#endif
#ifdef VMU
                Q->nbt[s]  += q->nbt[s];
                Q->nbt[s]  /= 2;
                Q->nbx[s]  += q->nbx[s];
                Q->nbx[s]  /= 2;
                Q->nby[s]  += q->nby[s];
                Q->nby[s]  /= 2;
                Q->nbn[s]  += q->nbn[s];
                Q->nbn[s]  /= 2;
#endif
            }
		}
	}
}

/**************************************************************************************************************************************************/

void regulateDissipativeCurrents(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ currrentVars, const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ u, int ncx, int ncy, int ncz) {
    
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);

#ifdef PIMUNU
				PRECISION pitt = currrentVars->pitt[s];
				PRECISION pitx = currrentVars->pitx[s];
				PRECISION pity = currrentVars->pity[s];
				PRECISION pitn = currrentVars->pitn[s];
				PRECISION pixx = currrentVars->pixx[s];
				PRECISION pixy = currrentVars->pixy[s];
				PRECISION pixn = currrentVars->pixn[s];
				PRECISION piyy = currrentVars->piyy[s];
				PRECISION piyn = currrentVars->piyn[s];
				PRECISION pinn = currrentVars->pinn[s];
#else
                PRECISION pitt = 0;
                PRECISION pitx = 0;
                PRECISION pity = 0;
                PRECISION pitn = 0;
                PRECISION pixx = 0;
                PRECISION pixy = 0;
                PRECISION pixn = 0;
                PRECISION piyy = 0;
                PRECISION piyn = 0;
                PRECISION pinn = 0;
#endif
#ifdef PI
				PRECISION Pi = currrentVars->Pi[s];
#else
				PRECISION Pi = 0;
#endif

				PRECISION ut = u->ut[s];
				PRECISION ux = u->ux[s];
				PRECISION uy = u->uy[s];
				PRECISION un = u->un[s];

				PRECISION xi0 = (PRECISION)(1.0);
				PRECISION rhomax = (PRECISION)(10.0);
				//PRECISION xi0 = (PRECISION)(0.1);
				//PRECISION rhomax = (PRECISION)(0.8);
				PRECISION t2 = t*t;
		//		PRECISION pipi = pitt*pitt-2*(pitx*pitx+pity*pity-pixy*pixy+t2*(pitn*pitn-pixn*pixn-piyn*piyn))+pixx*pixx+piyy*piyy+pinn*pinn*t2*t2;
				PRECISION pipi = pitt*pitt-2*pitx*pitx-2*pity*pity+pixx*pixx+2*pixy*pixy+piyy*piyy-2*pitn*pitn*t2+2*pixn*pixn*t2+2*piyn*piyn*t2+pinn*pinn*t2*t2;
		        if(isnan(pipi)==1) printf("found pipi Nan\n");
				PRECISION spipi = sqrt(fabs(pipi+3*Pi*Pi));
                if(isnan(spipi)==1) printf("found spipi Nan\n");
				PRECISION pimumu = pitt - pixx - piyy - pinn*t*t;
				PRECISION piu0 = -(pitn*t2*un) + pitt*ut - pitx*ux - pity*uy;
				PRECISION piu1 = -(pixn*t2*un) + pitx*ut - pixx*ux - pixy*uy;
				PRECISION piu2 = -(piyn*t2*un) + pity*ut - pixy*ux - piyy*uy;
				PRECISION piu3 = -(pinn*t2*un) + pitn*ut - pixn*ux - piyn*uy;

				PRECISION a1 = spipi/rhomax/sqrtf(e[s]*e[s]+3*p[s]*p[s]);
		        if(isnan(a1)==1) printf("found a1 Nan\n");
				PRECISION a2 = pimumu/xi0/rhomax/spipi;
				PRECISION a3 = piu0/xi0/rhomax/spipi;
				PRECISION a4 = piu1/xi0/rhomax/spipi;
				PRECISION a5 = piu2/xi0/rhomax/spipi;
				PRECISION a6 = piu3/xi0/rhomax/spipi;
				PRECISION a12 = fmax(a1,a2);
				PRECISION a34 = fmax(a3,a4);
				PRECISION a56 = fmax(a5,a6);
				PRECISION a3456 = fmax(a34,a56);
				PRECISION rho = fmax(a12,a3456);

				if(isnan(rho)==1) printf("found rho Nan\n");
				PRECISION fac = tanh(rho)/rho;
				if(fabs(rho)<1.e-7) fac = 1;
		//		PRECISION fac = 1;
				if(isnan(fac)==1) printf("found Nan\n");
		//		fac = 1;

#ifdef PIMUNU
				currrentVars->pitt[s] *= fac;
				currrentVars->pitx[s] *= fac;
				currrentVars->pity[s] *= fac;
				currrentVars->pitn[s] *= fac;
				currrentVars->pixx[s] *= fac;
				currrentVars->pixy[s] *= fac;
				currrentVars->pixn[s] *= fac;
				currrentVars->piyy[s] *= fac;
				currrentVars->piyn[s] *= fac;
				currrentVars->pinn[s] *= fac;
#endif
			}
		}
	}
}

/**************************************************************************************************************************************************/

void rungeKutta2(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q, void * latticeParams, void * hydroParams) {
    
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
    
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dy = (PRECISION)(lattice->latticeSpacingY);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

	PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);

	//===================================================
	// STEP 1:
	//===================================================
    
	eulerStepKernelSource(t, q, qS, e, ep, p, u, up, ncx, ncy, ncz, dt, dx, dy, dz, etabar, rhob, rhobp);
    eulerStepKernelX(t, q, qS, u, e, ncx, ncy, ncz, dt, dx, rhob);
	eulerStepKernelY(t, q, qS, u, e, ncx, ncy, ncz, dt, dy, rhob);
    eulerStepKernelZ(t, q, qS, u, e, ncx, ncy, ncz, dt, dz, rhob);
    
	t+=dt;

	setInferredVariablesKernel(qS, eS, p, uS, t, latticeParams, rhobS);

#ifndef IDEAL
	regulateDissipativeCurrents(t, qS, eS, p, uS, ncx, ncy, ncz);
#endif

	setGhostCells(qS, eS, p, uS, latticeParams, rhobS);

    
	//===================================================
	// STEP 2:
	//===================================================
    
	eulerStepKernelSource(t, qS, Q, eS, e, p, uS, u, ncx, ncy, ncz, dt, dx, dy, dz, etabar, rhobS, rhob);
	eulerStepKernelX(t, qS, Q, uS, eS, ncx, ncy, ncz, dt, dx, rhobS);
	eulerStepKernelY(t, qS, Q, uS, eS, ncx, ncy, ncz, dt, dy, rhobS);
	eulerStepKernelZ(t, qS, Q, uS, eS, ncx, ncy, ncz, dt, dz, rhobS);

	convexCombinationEulerStepKernel(q, Q, ncx, ncy, ncz);

	swapFluidVelocity(&up, &u);
    swapPrimaryVariables(&ep, &e);//Lipei
    swapPrimaryVariables(&rhobp, &rhob);//Lipei
    
	setInferredVariablesKernel(Q, e, p, u, t, latticeParams, rhob);
    
#ifndef IDEAL
	regulateDissipativeCurrents(t, Q, e, p, u, ncx, ncy, ncz);
#endif

	setGhostCells(Q, e, p, u, latticeParams, rhob);
}
