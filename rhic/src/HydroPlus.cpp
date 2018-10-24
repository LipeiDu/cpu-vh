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

/**************************************************************************************************************************************************/
/* Conductivity, heat capacity, correlation length and relaxation coefficients of different slow modes, and their derivatives
/**************************************************************************************************************************************************/

#define lambdaT 1.0

// heat capacity density
PRECISION Cp(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// correlation length
PRECISION xi(PRECISION e, PRECISION rhob){
    return 1.0; // correlation length, to be fixed
}

// derivative of heat capacity density with respect to energy
PRECISION dCpde(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// derivative of correlation length with respect to energy
PRECISION dxide(PRECISION e, PRECISION rhob){
    return 1.0; // correlation length, to be fixed
}

// derivative of heat capacity density with respect to baryon density
PRECISION dCpdrhob(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// derivative of correlation length with respect to baryon density
PRECISION dxidrhob(PRECISION e, PRECISION rhob){
    return 1.0; // correlation length, to be fixed
}

// universal function
PRECISION f2(PRECISION x){
    return 1.0 / (1.0 + x*x); // Eq. (93)
}

//relaxation coefficents of fluctuations
PRECISION relaxationCoefficientPhiQ(PRECISION e, PRECISION rhob, PRECISION Q)
{
    PRECISION corrL = xi(e, rhob);
    PRECISION qL = Q * corrL;
    PRECISION qL2 = qL * qL;
    PRECISION qL3 = qL2 * qL;
    PRECISION qL4 = qL3 * qL;
    
    return lambdaT/(Cp(e, rhob) * corrL * corrL) * (qL2 + qL4);
}


/**************************************************************************************************************************************************/
/* slow modes out and at equilibrium with different Q
/**************************************************************************************************************************************************/

// slow modes with Q=0
PRECISION equilibriumPhi0(PRECISION e, PRECISION rhob){
    return Cp(e, rhob) / (rhob * rhob); // Cp needs to be fixed, Eq.(90)
}

// slow modes with Q
PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION Q)
{
    PRECISION corrL = xi(e, rhob);
    PRECISION qL = Q * corrL;
    return equilibriumPhi0(e, rhob) * f2(qL); // Magnitude of mode Q at Equilibrium, Eq. (89)
}

// initialization of slow modes
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
    
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    PRECISION es = e[s];
                    PRECISION rhobs = rhob[s];
                    
                    eqPhiQ->phiQ[n][s] = equilibriumPhiQ(es, rhobs, Qvec[n]);
                    eqPhiQp->phiQ[n][s] = eqPhiQ->phiQ[n][s];
                    
                    q->phiQ[n+ALL_NUMBER_CONSERVED_VARIABLES][s] = eqPhiQ->phiQ[n][s];
                    Q->phiQ[n+ALL_NUMBER_CONSERVED_VARIABLES][s] = eqPhiQ->phiQ[n][s];
                }
            }
        }
    }
}

/**************************************************************************************************************************************************/
/* contributions from slow modes to inferred variables for a specific cell
/**************************************************************************************************************************************************/

// note: the integrands of alpha and beta only work for f2 defined above.
void getPrimaryVariablesFromSlowModes(PRECISION * const __restrict__ deltaS, PRECISION * const __restrict__ deltaAlpha, PRECISION * const __restrict__ deltaBeta, const PRECISION * const __restrict__ equiPhiQ, const PRECISION * const __restrict__ PhiQ, PRECISION e, PRECISION rhob)
{

    PRECISION entropy = 0.0;
    PRECISION alpha = 0.0;
    PRECISION beta = 0.0;
    
    PRECISION corrL = xi(e, rhob); // correlation length
    PRECISION heatC = Cp(e, rhob); // heat capacity
    PRECISION dCp_de = dCpde(e, rhob);
    PRECISION dCp_drhob = dCpdrhob(e, rhob);
    PRECISION dxi_de = dxide(e, rhob);
    PRECISION dxi_drhob = dxidrhob(e, rhob);
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        PRECISION deltaQ = Qvec[n+1] - Qvec[n]; // to be fixed
        PRECISION Q2 = Qvec[n]*Qvec[n];

        PRECISION qL = Qvec[n] * corrL;
        PRECISION qL2= qL * qL;
        
        PRECISION eqPhiQ = equiPhiQ[n];
        PRECISION phiRatio = PhiQ[n] / equiPhiQ[n];
        PRECISION deltaPhiQ = PhiQ[n] - equiPhiQ[n];
        PRECISION eqPhiQ2 = equiPhiQ[n] * equiPhiQ[n];
        PRECISION fac = deltaPhiQ / eqPhiQ2;
        PRECISION facC = eqPhiQ / heatC;
        PRECISION facX = 2 * Q2 * corrL * eqPhiQ2 * rhob * rhob / heatC;
        
        // delta beta, Eq.(106)
        PRECISION deqPhiQde =  facC * dCp_de - facX * dxi_de;
        PRECISION intBeta = Q2 * fac * deqPhiQde;
        beta += deltaQ * intBeta;
        
        // delta alpha, Eq.(106)
        PRECISION deqPhiQdrhob =  facC * dCp_drhob - facX * dxi_drhob - 2 / rhob * eqPhiQ;
        PRECISION intAlpha = Q2 * fac * deqPhiQdrhob;
        alpha += deltaQ * intAlpha;
        
        // delta entropy, Eq.(85)
        PRECISION intEntropy = Q2 * (log(phiRatio) - phiRatio + 1);
        entropy += deltaQ * intEntropy;
    }
    
    *deltaS = entropy;
    *deltaAlpha = alpha;
    *deltaBeta = beta;
}
