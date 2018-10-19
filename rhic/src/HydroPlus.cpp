//
//  HydroPlus.cpp
//  
//
//  Created by Lipei Du on 10/17/18.
//

// Equation indices from PRD 98 (2018) 036006

#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "../include/InitialConditions.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"
#include "../include/EquationOfState.h"
#include "../include/HydroPlus.h"

#define HydroPlus

#define lambdaT 1.0


PRECISION *Qvec; // Q vectors of slow modes
SLOW_MODES *PhiQ, *PhiQp, *PhiQS; // Slow modes: updated, previous, intermediate values
SLOW_MODES *eqPhiQ, *eqPhiQp, *eqPhiQS; // Slow modes at equilibrium: updated, previous, intermediate values

// allocate memory
void allocateHostMemorySlowModes(int len) {
    
    size_t bytes = sizeof(PRECISION);

    Qvec = (PRECISION *)calloc(NUMBER_SLOW_MODES, bytes);
    
    PhiQ  = (SLOW_MODES *)calloc(1, sizeof(SLOW_MODES));
    PhiQp = (SLOW_MODES *)calloc(1, sizeof(SLOW_MODES));
    PhiQS = (SLOW_MODES *)calloc(1, sizeof(SLOW_MODES));
    
    eqPhiQ  = (SLOW_MODES *)calloc(1, sizeof(SLOW_MODES));
    eqPhiQp = (SLOW_MODES *)calloc(1, sizeof(SLOW_MODES));
    eqPhiQS = (SLOW_MODES *)calloc(1, sizeof(SLOW_MODES));
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        PhiQ->phiQ[n]  = (PRECISION *)calloc(len, bytes);
        PhiQp->phiQ[n] = (PRECISION *)calloc(len, bytes);
        PhiQS->phiQ[n] = (PRECISION *)calloc(len, bytes);
        
        eqPhiQ->phiQ[n]  = (PRECISION *)calloc(len, bytes);
        eqPhiQp->phiQ[n] = (PRECISION *)calloc(len, bytes);
        eqPhiQS->phiQ[n] = (PRECISION *)calloc(len, bytes);
    }
}


PRECISION Cp(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// correlation length
PRECISION xi(PRECISION e, PRECISION rhob){
    return 1.0; // correlation length, to be fixed
}

// universal function
PRECISION f2(PRECISION x){
    return 1.0 / (1.0 + x*x); // Eq. (93)
}

// slow modes with Q=0
PRECISION equilibriumPhi0(PRECISION e, PRECISION rhob){
    return Cp(e, rhob) / (rhob * rhob); // Cp needs to be fixed, Eq.(90)
}

// slow modes with Q
PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION Q){
    PRECISION corrL = xi(e, rhob);
    PRECISION qL = Q * corrL;
    return equilibriumPhi0(e, rhob) * f2(qL); // Magnitude of mode Q at Equilibrium, Eq. (89)
}

//relaxation of fluctuations
PRECISION relaxationCoefficientPhiQ(PRECISION e, PRECISION rhob, PRECISION Q){
    PRECISION corrL = xi(e, rhob);
    PRECISION qL = Q * corrL;
    PRECISION qL2 = qL * qL;
    PRECISION qL3 = qL2 * qL;
    PRECISION qL4 = qL3 * qL;
    
    return lambdaT/(Cp(e, rhob) * corrL * corrL) * (qL2 + qL4);
}

void setInitialConditionSlowModes(void * latticeParams, void * hydroParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    Qvec[0] = 0.1;
    Qvec[1] = 1.0;
    Qvec[2] = 2.0;
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        PRECISION Qn = Qvec[n];
        
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    
                    PRECISION es = e[s];
                    PRECISION rhobs = rhob[s];
                
                    eqPhiQ->phiQ[n][s] = equilibriumPhiQ(es, rhobs, Qn);
                    eqPhiQp->phiQ[n][s] = eqPhiQ->phiQ[n][s];
                
                    PhiQ->phiQ[n][s] = eqPhiQ->phiQ[n][s];
                    PhiQp->phiQ[n][s] = eqPhiQ->phiQ[n][s];
                }
            }
        }
    }
}


//contributions from slow modes to primary variables
PRECISION entropySlowModes(const PRECISION * const __restrict__ equilibriumPhiQ, const PRECISION * const __restrict__ PhiQ, const PRECISION * const __restrict__ Q){
    // entropy of slow modes on a single cell
    
    PRECISION entropy = 0.0;
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        PRECISION deltaQ = Q[n+1] - Q[n]; // to be fixed
        PRECISION phiRatio = PhiQ[n] / equilibriumPhiQ[n];
        PRECISION integrand = Q[n]*Q[n] * (log(phiRatio) - phiRatio + 1);
        entropy += deltaQ * integrand;
    }
    
    return entropy;
}


//source terms of extra modes, Eq. (76)
void setSlowModesSourceTerms(PRECISION * const __restrict__ PhiQRHS, const PRECISION * const __restrict__ equilibriumPhiQ, const PRECISION * const __restrict__ PhiQ, const PRECISION * const __restrict__ Qvec, PRECISION e, PRECISION rhob, PRECISION ut, PRECISION dkvk){
    
    PRECISION gammaQ[NUMBER_SLOW_MODES];
    PRECISION utInv = 1.0/ut;
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        PRECISION Qn = Qvec[n];
        gammaQ[n] = relaxationCoefficientPhiQ(e, rhob, Qn);
        PhiQRHS[n] = utInv * (-2) * gammaQ[n] * (PhiQ[n] - equilibriumPhiQ[n]) + PhiQ[n] * dkvk;
    }
    
}

void loadSlowModesSourceTerms(const SLOW_MODES * const __restrict__ equilibriumSlowModes, const SLOW_MODES * const __restrict__ SlowModes, PRECISION * const __restrict__ S, const PRECISION * const __restrict__ Qvec, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ rhobvec, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz){
    // Modes, slow modes; S, source terms;
    
    //=========================================================
    // primary variables
    //=========================================================
    
    PRECISION e = evec[s];
    PRECISION rhob = rhobvec[s];
    
    PRECISION *utvec = u->ut;
    PRECISION *uxvec = u->ux;
    PRECISION *uyvec = u->uy;
    PRECISION *unvec = u->un;
    
    PRECISION ut = utvec[s];
    PRECISION ux = uxvec[s];
    PRECISION uy = uyvec[s];
    PRECISION un = unvec[s];
    
    PRECISION *PhiQvec[NUMBER_SLOW_MODES];
    PRECISION PhiQ[NUMBER_SLOW_MODES];
    PRECISION *equilibriumPhiQvec[NUMBER_SLOW_MODES];
    PRECISION equilibriumPhiQ[NUMBER_SLOW_MODES];
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        equilibriumPhiQvec[n] = equilibriumSlowModes->phiQ[n];
        equilibriumPhiQ[n] = equilibriumPhiQvec[n][s];
        
        PhiQvec[n] = SlowModes->phiQ[n];
        PhiQ[n] = PhiQvec[n][s];
        
    }
        
    //=========================================================
    // spatial derivatives of primary variables
    //=========================================================
    
    PRECISION facX = 1/d_dx/2;
    PRECISION facY = 1/d_dy/2;
    PRECISION facZ = 1/d_dz/2;
    
    PRECISION dxut = (*(utvec + s + 1) - *(utvec + s - 1)) * facX;
    PRECISION dxux = (*(uxvec + s + 1) - *(uxvec + s - 1)) * facX;
    
    PRECISION dyut = (*(utvec + s + d_ncx) - *(utvec + s - d_ncx)) * facY;
    PRECISION dyuy = (*(uyvec + s + d_ncx) - *(uyvec + s - d_ncx)) * facY;
    
    int stride = d_ncx * d_ncy;
    PRECISION dnut = (*(utvec + s + stride) - *(utvec + s - stride)) * facZ;
    PRECISION dnun = (*(unvec + s + stride) - *(unvec + s - stride)) * facZ;
    
    PRECISION vx = ux/ut;
    PRECISION vy = uy/ut;
    PRECISION vn = un/ut;
    PRECISION dxvx = (dxux - vx * dxut)/ ut;
    PRECISION dyvy = (dyuy - vy * dyut)/ ut;
    PRECISION dnvn = (dnun - vn * dnut)/ ut;
    PRECISION dkvk = dxvx + dyvy + dnvn;
    
    //=========================================================
    // extra modes source terms
    //=========================================================
    
    PRECISION PhiQRHS[NUMBER_SLOW_MODES];
    
    setSlowModesSourceTerms(PhiQRHS, equilibriumPhiQ, PhiQ, Qvec, e, rhob, ut, dkvk);
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n) S[n] = PhiQRHS[n];
    
}
