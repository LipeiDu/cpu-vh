/*
 * DynamicalVariables.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_


#define NUMBER_CONSERVATION_LAWS 4

#define PIMUNU
#define PI

//Lipei
#define NBMU
#define VMU

#ifndef NBMU
#define NUMBER_BARYON_COMPONENTS 0
#else
#define NUMBER_BARYON_COMPONENTS 1
#endif

#define ALL_NUMBER_CONSERVATION_LAWS (NUMBER_CONSERVATION_LAWS+NUMBER_BARYON_COMPONENTS)
//Lipei

/*********************************************************/
#ifndef PI
#define NUMBER_PI_COMPONENTS 0
#else
#define NUMBER_PI_COMPONENTS 1
#endif
// bulk pressure; Lipei's comment

#ifndef VMU
#define NUMBER_PROPAGATED_VMU_COMPONENTS 0
#else
#define NUMBER_PROPAGATED_VMU_COMPONENTS 4
#endif
// baryon diffusion current; Lipei's comment

#ifndef PIMUNU
#define NUMBER_PROPAGATED_PIMUNU_COMPONENTS 0
#else
#define NUMBER_PROPAGATED_PIMUNU_COMPONENTS 10
#endif
// shear stress tensor; Lipei's comment

#define NUMBER_DISSIPATIVE_CURRENTS (NUMBER_PI_COMPONENTS+NUMBER_PROPAGATED_PIMUNU_COMPONENTS)

#define ALL_NUMBER_DISSIPATIVE_CURRENTS (NUMBER_PI_COMPONENTS+NUMBER_PROPAGATED_VMU_COMPONENTS+NUMBER_PROPAGATED_PIMUNU_COMPONENTS)
//Lipei

#if ALL_NUMBER_DISSIPATIVE_CURRENTS==0
#define IDEAL
#endif

#define NUMBER_CONSERVED_VARIABLES (NUMBER_CONSERVATION_LAWS+NUMBER_DISSIPATIVE_CURRENTS)
//excluding baryon part; Lipei's comment

#define ALL_NUMBER_CONSERVED_VARIABLES (ALL_NUMBER_CONSERVATION_LAWS+ALL_NUMBER_DISSIPATIVE_CURRENTS)
//by Lipei; baryon part is defined seperately so that arries can be taken care of more easily.

/*********************************************************/

#define PRECISION double

typedef struct
{
	PRECISION *ttt;
	PRECISION *ttx;
	PRECISION *tty;
	PRECISION *ttn;
#ifdef PIMUNU
	PRECISION *pitt;
	PRECISION *pitx;
	PRECISION *pity;
	PRECISION *pitn;
	PRECISION *pixx;
	PRECISION *pixy;
	PRECISION *pixn;
	PRECISION *piyy;
	PRECISION *piyn;
	PRECISION *pinn;
#endif
#ifdef PI
	PRECISION *Pi;
#endif
//baryon; Lipei
#ifdef NBMU
    PRECISION *Nbt;
#endif
#ifdef VMU
    PRECISION *nbt;
    PRECISION *nbx;
    PRECISION *nby;
    PRECISION *nbn;
#endif
} CONSERVED_VARIABLES;

typedef struct
{
	PRECISION *ut;
	PRECISION *ux;
	PRECISION *uy;
	PRECISION *un;
} FLUID_VELOCITY;

typedef struct
{
    PRECISION *sourcet;
    PRECISION *sourcex;
    PRECISION *sourcey;
    PRECISION *sourcen;
    PRECISION *sourceb;
} DYNAMICAL_SOURCE;//lipei

extern CONSERVED_VARIABLES *q,*Q,*qS;
extern FLUID_VELOCITY *u,*up,*uS,*uSS;
extern PRECISION *e, *p;
extern PRECISION *rhob;//Lipei
extern DYNAMICAL_SOURCE *Source;//Lipei

int columnMajorLinearIndex(int i, int j, int k, int nx, int ny);

void allocateHostMemory(int len);

void setConservedVariables(double t, void * latticeParams);
void setCurrentConservedVariables();
void swapFluidVelocity(FLUID_VELOCITY **arr1, FLUID_VELOCITY **arr2) ;

void setGhostCells(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob//by Lipei
);

void setGhostCellsKernelI(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob//by Lipei
);

void setGhostCellsKernelJ(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob//by Lipei
);

void setGhostCellsKernelK(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob//by Lipei
);

void freeHostMemory();

#endif /* DYNAMICALVARIABLES_H_ */
