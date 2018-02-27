/*
 * InitialConditions.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <math.h> // for math functions
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP

#include <iostream>// by Lipei
#include <fstream>// by Lipei
#include <iomanip>//by Lipei

#include "edu/osu/rhic/trunk/ic/InitialConditions.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/core/ic/GlauberModel.h"
#include "edu/osu/rhic/core/ic/MonteCarloGlauberModel.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/trunk/eos/EquationOfState.h"

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)

using namespace std;//Lipei


//*********************************************************************************************************\
//* Read in all initial profiles from a single or seperate file
//*********************************************************************************************************/

//this reads all hydro variables from a single file; this way we do not need to fetch the coordinates many times
//note that the file must contain values for all dissipative currents, even if they are zero !!!
void setInitialTmunuFromFile(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    float x, y, z, e_in, p_in, ut_in, ux_in, uy_in, un_in;
    //#ifdef PIMUNU
    float pitt_in, pitx_in, pity_in, pitn_in, pixx_in, pixy_in, pixn_in, piyy_in, piyn_in, pinn_in;
    //#endif
    //#ifdef PI
    float Pi_in;
    //#endif
    FILE *fileIn;
    char fname[255];
    
    sprintf(fname, "%s/%s", rootDirectory, "/input/Tmunu.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open Tmunu.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &x, &y, &z, &e_in, &p_in, &ut_in, &ux_in, &uy_in, &un_in, &pitt_in, &pitx_in, &pity_in, &pitn_in, &pixx_in, &pixy_in, &pixn_in, &piyy_in, &piyn_in, &pinn_in, &Pi_in);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    e[s] =  (PRECISION) e_in;
                    //ep[s] = (PRECISION) e_in; //set previous step to same value
                    p[s] = p_in;
                    u->ut[s] = ut_in;
                    u->ux[s] = ux_in;
                    u->uy[s] = uy_in;
                    u->un[s] = un_in;
                    up->ut[s] = ut_in; //set previous step to same value
                    up->ux[s] = ux_in; //...
                    up->uy[s] = uy_in;
                    up->un[s] = un_in;
#ifdef PIMUNU
                    q->pitt[s] = pitt_in;
                    q->pitx[s] = pitx_in;
                    q->pity[s] = pity_in;
                    q->pitn[s] = pitn_in;
                    q->pixx[s] = pixx_in;
                    q->pixy[s] = pixy_in;
                    q->pixn[s] = pixn_in;
                    q->piyy[s] = piyy_in;
                    q->piyn[s] = piyn_in;
                    q->pinn[s] = pinn_in;
#endif
#ifdef PI
                    q->Pi[s] = Pi_in;
#endif
                }
            }
        }
    }
    fclose(fileIn);
}

//this function reads a separate file for every hydrodynamic variable
void setInitialTmunuFromFiles(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    float x, y, z, value;
    FILE *fileIn;
    char fname[255];
    
    //energy density
    sprintf(fname, "%s/%s", rootDirectory, "/input/e.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open e.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    e[s] =  (PRECISION) value;
                    //ep[s] = (PRECISION) value;
                    //printf("e [ %d ] = %f\n", s, e[s]);
                }
            }
        }
    }
    fclose(fileIn);
    
    //pressure
    sprintf(fname, "%s/%s", rootDirectory, "/input/p.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open p.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    p[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //ut
    sprintf(fname, "%s/%s", rootDirectory, "/input/ut.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ut.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->ut[s] =  (PRECISION) value;
                    up->ut[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //ux
    sprintf(fname, "%s/%s", rootDirectory, "/input/ux.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ux.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->ux[s] =  (PRECISION) value;
                    up->ux[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //uy
    sprintf(fname, "%s/%s", rootDirectory, "/input/uy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open uy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->uy[s] =  (PRECISION) value;
                    up->uy[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //un
    sprintf(fname, "%s/%s", rootDirectory, "/input/un.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open un.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->un[s] =  (PRECISION) value;
                    up->un[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
#ifdef PIMUNU
    //pitt
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitt.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitt.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitt[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pitx
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitx[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pity
    sprintf(fname, "%s/%s", rootDirectory, "/input/pity.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pity.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pity[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pitn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pixx
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixx[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pixy
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixy[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pixn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //piyy
    sprintf(fname, "%s/%s", rootDirectory, "/input/piyy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open piyy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->piyy[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //piyn
    sprintf(fname, "%s/%s", rootDirectory, "/input/piyn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open piyn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->piyn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pinn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pinn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pinn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pinn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
#ifdef PI
    //bulk
    sprintf(fname, "%s/%s", rootDirectory, "/input/bulk.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open bulk.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->Pi[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
}

//*********************************************************************************************************\
//* Set initial flow profile
//*	- u^\mu = (1, 0, 0, 0)
//* 	- No transverse flow (ux = uy = 0)
//*	- Longitudinal scaling flow (u_z = z/t, i.e. un = 0)
//*********************************************************************************************************/
void setFluidVelocityInitialCondition(void * latticeParams, void * hydroParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double t0 = hydro->initialProperTimePoint;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				PRECISION ux = 0;
				PRECISION uy = 0;
				PRECISION un = 0;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = sqrt(1+ux*ux+uy*uy+t0*t0*un*un);
                
                //intialize the flow velocity of the previous step; Lipei
                up->ux[s] = 0;
                up->uy[s] = 0;
                up->un[s] = 0;
                up->ut[s] = sqrt(1+ux*ux+uy*uy+t0*t0*un*un);
			}
		}
	}
}

//*********************************************************************************************************\
//* Set initial baryon diffusion current//Lipei
//*********************************************************************************************************/
void setbnmuInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    printf("Initialize \\nb^\\mu to be zero.\n");
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
#ifdef VMU
                q->nbt[s] = 0;
                q->nbx[s] = 0;
                q->nby[s] = 0;
                q->nbn[s] = 0;
#endif
            }
        }
    }
    return;
}

//*********************************************************************************************************\
//* Set initial shear-stress tensor \pi^\mu\nu
//*	- Navier-Stokes value, i.e. \pi^\mu\nu = 2 * (\epsilon + P) / T * \eta/S * \sigma^\mu\nu
//* 	- No initial pressure anisotropies (\pi^\mu\nu = 0)
//*********************************************************************************************************/
void setPimunuNavierStokesInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

	PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);
	PRECISION t = hydro->initialProperTimePoint;

	PRECISION e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				PRECISION T = effectiveTemperature(e[s], rhob[s]);
				if (T == 0) T = 1.e-3;
				PRECISION pinn = -4.0/(3*t*t*t)*etabar*(e[s]+p[s])/T;//was wrong by factor of 2
#ifdef PIMUNU
				q->pitt[s] = 0;
				q->pitx[s] = 0;
				q->pity[s] = 0;
				q->pitn[s] = 0;
				q->pixx[s] = -t*t*pinn/2;
				q->pixy[s] = 0;
				q->pixn[s] = 0;
				q->piyy[s] = -t*t*pinn/2;
				q->piyn[s] = 0;
				q->pinn[s] = pinn;
#endif
#ifdef PI
#define A_1 -13.77
#define A_2 27.55
#define A_3 13.45

#define LAMBDA_1 0.9
#define LAMBDA_2 0.25
#define LAMBDA_3 0.9
#define LAMBDA_4 0.22

#define SIGMA_1 0.025
#define SIGMA_2 0.13
#define SIGMA_3 0.0025
#define SIGMA_4 0.022
				PRECISION x = T/1.01355;
				PRECISION zetabar = A_1*x*x + A_2*x - A_3;
				if(x > 1.05)
					zetabar = LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
				else if(x < 0.995)
					zetabar = LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
				q->Pi[s] = -zetabar*(e[s]+p[s])/T/t;
#endif
			}
		}
	}
}

void setPimunuInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int initializePimunuNavierStokes = hydro->initializePimunuNavierStokes;
	if (initializePimunuNavierStokes==1) {
		printf("Initialize \\pi^\\mu\\nu to its asymptotic Navier-Stokes value.\n");
#ifdef PI
		printf("Initialize \\Pi to its asymptotic Navier-Stokes value.\n");
#endif
		setPimunuNavierStokesInitialCondition(latticeParams, initCondParams, hydroParams);
		return;
	}
	else {
		printf("Initialize \\pi^\\mu\\nu to zero.\n");
		struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
		int nx = lattice->numLatticePointsX;
		int ny = lattice->numLatticePointsY;
		int nz = lattice->numLatticePointsRapidity;
		for(int i = 2; i < nx+2; ++i) {
			for(int j = 2; j < ny+2; ++j) {
				for(int k = 2; k < nz+2; ++k) {
					int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
#ifdef PIMUNU
			  		q->pitt[s] = 0;
			  		q->pitx[s] = 0;
			  		q->pity[s] = 0;
			  		q->pitn[s] = 0;
			  		q->pixx[s] = 0;
			  		q->pixy[s] = 0;
			  		q->pixn[s] = 0;
			  		q->piyy[s] = 0;
			  		q->piyn[s] = 0;
			  		q->pinn[s] = 0;
#endif
#ifdef PI
			  		q->Pi[s] = 0;
#endif
				}
			}
		}
		return;
	}
}

//*********************************************************************************************************\
//* Longitudinal initial energy density distribution
//*********************************************************************************************************/
void longitudinalEnergyDensityDistribution(double * const __restrict__ eL, void * latticeParams, void * initCondParams) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nz = lattice->numLatticePointsRapidity;
    
    double dz = lattice->latticeSpacingRapidity;
    
    double etaFlat = initCond->rapidityMean;
    double etaVariance = initCond->rapidityVariance;
    
    for(int k = 0; k < nz; ++k) {
        double eta = (k - (nz-1)/2)*dz;
        double etaScaled = fabs(eta) - etaFlat/2;
        double arg = -etaScaled * etaScaled / etaVariance / 2 * THETA_FUNCTION(etaScaled);
        eL[k] = exp(arg);
    }
}

//*********************************************************************************************************\
//* Longitudinal initial baryon density distribution
//*********************************************************************************************************/
void longitudinalBaryonDensityDistribution(double * const __restrict__ rhoLa, double * const __restrict__ rhoLb, void * latticeParams, void * initCondParams) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nz = lattice->numLatticePointsRapidity;
    double dz = lattice->latticeSpacingRapidity;
    
    double etaVariance1 = 0.1;
    double etaVariance2 = 1.0;
    double etaMean = 2.0;
    
    for(int k = 0; k < nz; ++k) {
        double eta = (k - (nz-1)/2)*dz;
        rhoLa[k] = 1/sqrt(2*3.1415926)*(1/etaVariance1*exp(-(eta-etaMean)*(eta-etaMean)/(2*etaVariance1*etaVariance1))*THETA_FUNCTION(eta-etaMean)
                                     + 1/etaVariance2*exp(-(eta-etaMean)*(eta-etaMean)/(2*etaVariance2*etaVariance2))*THETA_FUNCTION(etaMean-eta));
        rhoLb[k] = 1/sqrt(2*3.1415926)*(1/etaVariance1*exp(-(-eta-etaMean)*(-eta-etaMean)/(2*etaVariance1*etaVariance1))*THETA_FUNCTION(-eta-etaMean)
                                     + 1/etaVariance2*exp(-(-eta-etaMean)*(-eta-etaMean)/(2*etaVariance2*etaVariance2))*THETA_FUNCTION(etaMean+eta));
    }
}


//*********************************************************************************************************\
//* Constant initial energy density distribution
//*********************************************************************************************************/
void setConstantEnergyDensityInitialCondition(void * latticeParams, void * initCondParams) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	double initialEnergyDensity = initCond->initialEnergyDensity;
    double initialBaryonDensity = initCond->initialBaryonDensity;//Lipei

	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
    
    double T0 = 3.05;
    double ed = equilibriumEnergyDensity(T0);
    
    double rhobd = initialBaryonDensity;//Lipei
    double rhoLa[nz], rhoLb[nz];//Lipei
    longitudinalBaryonDensityDistribution(rhoLa, rhoLb, latticeParams, initCondParams);//Lipei

    char rhobb[] = "output/baryon_density.dat";//Lipei
    ofstream baryondens(rhobb);//Lipei
    
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				e[s] = (PRECISION) ed;
                rhob[s] = rhoLa[k-2] * rhobd + rhoLb[k-2] * rhobd + 1.e-4; //Lipei
                p[s] = equilibriumPressure(e[s], rhob[s]);
                
                //ep[s] = e[s];
                //rhobp[s] = rhob[s];
                
                baryondens << setprecision(3) << setw(5) << i <<setprecision(3)  << setw(5)  << j
                << setprecision(3) << setw(5) << k << setprecision(6) << setw(18) << rhob[s] << endl;//Lipei
			}
		}
	}
    
    baryondens.close();//Lipei
    printf("Baryon density is initialized.\n");
}


//*********************************************************************************************************\
//* Continuous optical glauber Glauber initial energy density distribution
//*********************************************************************************************************/
void setGlauberInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
	double T0 = 3.05;
	e0 = (double) equilibriumEnergyDensity(T0);

	double eT[nx*ny], eL[nz];
    double rhoLa[nz], rhoLb[nz];//Lipei
    double Ta[nx*ny], Tb[nx*ny];//Lipei
    
    longitudinalBaryonDensityDistribution(rhoLa, rhoLb, latticeParams, initCondParams);//Lipei
	energyDensityTransverseProfileAA(eT, nx, ny, dx, dy, initCondParams, Ta, Tb);
	longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);

    char rhobb[] = "output/baryon_density.dat";//Lipei
    ofstream baryondens(rhobb);//Lipei
    
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			double energyDensityTransverse = e0 * eT[i-2+(j-2)*nx];
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				double energyDensityLongitudinal = eL[k-2];
                double ed = (energyDensityTransverse * energyDensityLongitudinal);// + 1.e-3;
				e[s] = (PRECISION) ed;
                rhob[s] = rhoLa[k-2]*Ta[i-2+(j-2)*nx] + rhoLb[k-2]*Tb[i-2+(j-2)*nx];// + 1.e-4; //Lipei
                p[s] = equilibriumPressure(e[s], rhob[s]);
                
                //ep[s] = e[s];
                //rhobp[s] = rhob[s];
                
                baryondens << setprecision(3) << setw(5) << i <<setprecision(3)  << setw(5)  << j
                           << setprecision(3) << setw(5) << k << setprecision(6) << setw(18) << rhob[s] << endl;//Lipei
			}
		}
	}

    baryondens.close();//Lipei
    printf("Baryon density is initialized.\n");
}

//*********************************************************************************************************\
//* Monte carlo Glauber initial energy density distribution
//*********************************************************************************************************/
void setMCGlauberInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
	double T0 = 2.03;
	e0 = (double) equilibriumEnergyDensity(T0);
    
    double t0 = hydro->initialProperTimePoint;//Lipei
    double Prefactor = 1/t0;//Lipei

    double eT[nx*ny], eL[nz];
    double rhoLa[nz], rhoLb[nz];//Lipei
    double Ta[nx*ny], Tb[nx*ny];//Lipei
    int n1, n2;//participants, Lipei
    
    monteCarloGlauberEnergyDensityTransverseProfile(eT, nx, ny, dx, dy, initCondParams, Ta, Tb, &n1, &n2);
    longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);
    longitudinalBaryonDensityDistribution(rhoLa, rhoLb, latticeParams, initCondParams);//Lipei
    
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				double energyDensityLongitudinal = eL[k-2];
                double eta = (k-2 - (nz-1)/2)*dz;//Lipei
                e[s] = 1.5*e0 * Prefactor * energyDensityLongitudinal * (Ta[i-2 + (j-2)*nx]*(2.0+eta)/4.0 + Tb[i-2 + (j-2)*nx]*(2.0-eta)/4.0) + 0.004;
                rhob[s] = 0.5*Prefactor * (n1 * rhoLa[k-2] * Ta[i-2+(j-2)*nx] + n2 * rhoLb[k-2] * Tb[i-2+(j-2)*nx]) + 0.0001; //Lipei
                p[s] = equilibriumPressure(e[s], rhob[s]);
			}
		}
	}
    printf("Baryon density is initialized.\n");
}

//*********************************************************************************************************\
//* Initial conditions for the Gubser ideal hydro test
//*		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
//*********************************************************************************************************/
void setIdealGubserInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			double T = 1.9048812623618392/pow(1 + pow(1 - pow(x,2) - pow(y,2),2) + 2*(1 + pow(x,2) + pow(y,2)),0.3333333333333333);
			double r = sqrt(x*x+y*y);
			double phi = atanh(2*1*r/(1+1+x*x+y*y));

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				e[s] = (PRECISION) (e0 * pow(T,4));
				p[s] = e[s]/3;
				u->ux[s] = (PRECISION) (sinh(phi)*x/r);
				u->uy[s] = (PRECISION) (sinh(phi)*y/r);
				u->un[s] = 0;
				u->ut[s] = sqrt(1 + u->ux[s]*u->ux[s] + u->uy[s]*u->uy[s]);
			}
		}
	}
}

//*********************************************************************************************************\
//* Initial conditions for the Gubser viscous hydro test
//*		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
//*********************************************************************************************************/
void setISGubserInitialCondition(void * latticeParams, const char *rootDirectory) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double x,y,ed,u1,u2,pitt,pitx,pity,pixx,pixy,piyy,pinn;

	FILE *file;
	char fname[255];
	sprintf(fname, "%s/%s", rootDirectory, "/rhic/rhic-trunk/src/test/resources/gubser/viscous/gubserIC.dat");
	file = fopen(fname, "r");

	double pitn=0;
	double pixn=0;
	double piyn=0;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			int status = fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
		    		&x,&y,&ed,&u1,&u2,&pixx,&piyy,&pixy,&pitt,&pitx,&pity,&pinn);
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				e[s] = (PRECISION) ed;
				p[s] = e[s]/3;
				u->ux[s] = u1;
				u->uy[s] = u2;
				u->un[s] = 0;
				u->ut[s] = sqrt(1 + u1*u1 + u2*u2);
#ifdef PIMUNU
        		q->pitt[s] = (PRECISION) pitt;
        		q->pitx[s] = (PRECISION) pitx;
        		q->pity[s] = (PRECISION) pity;
        		q->pitn[s] = (PRECISION) pitn;
        		q->pixx[s] = (PRECISION) pixx;
        		q->pixy[s] = (PRECISION) pixy;
        		q->pixn[s] = (PRECISION) pixn;
        		q->piyy[s] = (PRECISION) piyy;
        		q->piyn[s] = (PRECISION) piyn;
        		q->pinn[s] = (PRECISION) pinn;
#endif
			}
		}
	}
}

//*********************************************************************************************************\
//* Initial conditions for the relativistic Sod shock-tube test
//*		- set energy density, pressure, fluid velocity u^\mu
//*********************************************************************************************************/
void setSodShockTubeInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				if(x > 0) 	e[s] = (PRECISION) (0.00778147);
				else 			e[s] = (PRECISION) (0.124503);
//				if(y > 0) 	e[s] = (PRECISION) (0.00778147);
//				else 			e[s] = (PRECISION) (0.124503);
//				if(x > 0) 	e[s] = (PRECISION) (1.0);
//				else 			e[s] = (PRECISION) (100.0);
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

//*********************************************************************************************************\
//* Initial conditions for the relativistic Sod shock-tube test
//*		- set energy density, pressure, fluid velocity u^\mu
//*********************************************************************************************************/
void set2dSodShockTubeInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				if(y > x) 	e[s] = (PRECISION) (0.00778147);
//				if(atan(y/x)>0.7853981634) 	e[s] = (PRECISION) (0.00778147);
				else 			e[s] = (PRECISION) (0.124503);
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

//*********************************************************************************************************\
//* Initial conditions for the implosion in a box test
//*		- set energy density, pressure, fluid velocity u^\mu
//*********************************************************************************************************/
void setImplosionBoxInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
    double initialBaryonDensity = initCond->initialBaryonDensity;//Lipei

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
//				e[s] = (PRECISION) (0.00778147);
//				if (sqrt(x*x+y*y)<=0.15) e[s] = (PRECISION) (0.124503);

//				e[s] = (PRECISION) (0.124503);
//				if (sqrt(x*x+y*y)<=0.15) e[s] = (PRECISION) (0.00778147);
				e[s] = (PRECISION) (0.00778147);
				if (sqrt(x*x+y*y)<=0.15) e[s] = (PRECISION) (0.124503);

//				e[s] = (PRECISION) (1.0);
/*
				if (x < 1) {
					if (y < (1-x))
						e[s] = (PRECISION) (0.00778147);
				}
				else e[s] = (PRECISION) (0.124503);
//*/
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;

                rhob[s] = (PRECISION) initialBaryonDensity;//Lipei
			}
		}
	}
}

//*********************************************************************************************************\
//* Initial conditions for the implosion in a box test
//*		- set energy density, pressure, fluid velocity u^\mu
//*********************************************************************************************************/
void setRayleighTaylorInstibilityInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double Lx = ( (nx-1)/2.)*dx;
	double Ly = ( (ny-1)/2.)*dy;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				double gasGamma = 1.4;
				double gravity = 0.1;

				double rhoTop = 2.0;
				double rhoBot = 1.0;
				double pr0 = 0.01;
				double pert = 0.01;

				double yloc = 0.5 + pert*cos(M_PI*x);
				double pr;
				if(y > yloc) 	pr = rhoTop*gravity*(1-y);
				else 				pr = rhoTop*gravity*(1-yloc) + rhoBot*gravity*(yloc-y);
				e[s] = pr;
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

//*********************************************************************************************************\
//* Initial conditions for the implosion in a box test
//*		- set energy density, pressure, fluid velocity u^\mu
//*********************************************************************************************************/
void setGaussianPulseInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double Lx = ( (nx-1)/2.)*dx;
	double Ly = ( (ny-1)/2.)*dy;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				double xc = 0.5;
				double yc = 0.5;
			    double beta = 50.0;
			    double pr = 1.0 + 1e-1*exp(-beta*((x-xc)*(x-xc)+(y-yc)*(y-yc)));

				e[s] = (PRECISION) (pr/(1.4-1.));

				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = (1+cos(2*M_PI*x/Lx))*(1+cos(2*M_PI*y/Ly));
//				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = sqrt(1+u->uy[s]*u->uy[s]);
			}
		}
	}
}


//*********************************************************************************************************\
//* Initial conditions for the sound propagation test
//*        - set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
 /*********************************************************************************************************\

void setSoundPropagationInitialCondition(void * latticeParams, void * initCondParams) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    double initialEnergyDensity = initCond->initialEnergyDensity;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double T0 = 3.05; // -> 0.6 GeV
    double e0 = initialEnergyDensity*pow(T0, 4);
    
    double cs = 0.57735;
    double de = e0/100.;
    double PI = 3.141592653589793;
    double p0=e0/3;
    double lambda = (nx-1)*dx/2.;
    
    for(int i = 2; i < nx+2; ++i) {
        double x = (i-2 - (nx-1)/2.)*dx;
        double vx = cs*de/(e0+p0)*sin(2*PI*x/lambda);
        double ed = e0 + de*sin(2*PI*x/lambda);
        // periodic boundary conditions
        if (i==2) {
            vx = cs*de/(e0+p0)*sin(2*PI*abs(x)/lambda);
            ed = e0 + de*sin(2*PI*abs(x)/lambda);
        }
        
        double u0 = 1/sqrt(1-vx*vx);
        
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                e[s] = ed;
                p[s] = Pressure(e[s]);
                u->ux[s] = u0*vx;
                u->uy[s] = 0;
                u->un[s] = 0;
                u->ut[s] = u0;
                // initialize \pi^\mu\nu to zero
                q->pitt[s] = 0;
                q->pitx[s] = 0;
                q->pity[s] = 0;
                q->pitn[s] = 0;
                q->pixx[s] = 0;
                q->pixy[s] = 0;
                q->pixn[s] = 0;
                q->piyy[s] = 0;
                q->piyn[s] = 0;
                q->pinn[s] = 0;
            }
        }
    }
}*/

//*********************************************************************************************************\
//* Initial conditions to use.
//*	Set the energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny.
//* 	0 - constant energy density
//*		1 - Isreal-Stewart hydrodynamic Gubser flow test
//*		2 - Continous optical Glauber
//*		3 - Ideal hydrodynamic Gubser flow test
//*		4 - Monte carlo Glauber
//*		5 - Relativistic Sod shock-tube test
//*********************************************************************************************************/
void setInitialConditions(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int initialConditionType = initCond->initialConditionType;
	printf("Setting initial conditions: ");
	switch (initialConditionType) {
		case 0: {
			printf("Constant energy density.\n");
			setConstantEnergyDensityInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
            setbnmuInitialCondition(latticeParams, initCondParams, hydroParams);//Lipei
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 1: {
			printf("Isreal-Stewart hydrodynamic Gubser flow test.\n");
			setISGubserInitialCondition(latticeParams, rootDirectory);
			return;
		}
		case 2: {
			printf("Continous optical Glauber.\n");
			setGlauberInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
            setbnmuInitialCondition(latticeParams, initCondParams, hydroParams);//Lipei
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 3: {
			printf("Ideal hydrodynamic Gubser flow test.\n");
			setIdealGubserInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 4: {
			printf("Monte carlo Glauber.\n");
			setMCGlauberInitialCondition(latticeParams, initCondParams, hydroParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
            setbnmuInitialCondition(latticeParams, initCondParams, hydroParams);//Lipei
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 5: {
			printf("Relativistic Sod shock-tube test.\n");
			setSodShockTubeInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 6: {
			printf("Implosion in a box test.\n");
			setImplosionBoxInitialCondition(latticeParams, initCondParams);
            setbnmuInitialCondition(latticeParams, initCondParams, hydroParams);//Lipei
			return;
		}

		case 7: {
			printf("Rayleigh-Taylor instability test.\n");
			setRayleighTaylorInstibilityInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 8: {
			printf("Implosion in a box test.\n");
			setGaussianPulseInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 9: {
			printf("Relativistic 2d Sod shock-tube test.\n");
			set2dSodShockTubeInitialCondition(latticeParams, initCondParams);
			return;
		}
        case 10: {
            printf("Reading initial T ^mu nu from input/e.dat , input/p.dat , etc... \n");
            setInitialTmunuFromFiles(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
        }
        case 11: {
            printf("Reading initial T ^mu nu from /input/Tmunu.dat \n");
            setInitialTmunuFromFile(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
        }
		default: {
			printf("Initial condition type not defined. Exiting ...\n");
			exit(-1);
		}
	}
}
