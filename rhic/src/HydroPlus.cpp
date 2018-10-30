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

#define dQvec 0.2 // difference between Q vectors of slow modes

/**************************************************************************************************************************************************/
/* Conductivity, heat capacity, correlation length and relaxation coefficients of different slow modes, and their derivatives
/**************************************************************************************************************************************************/

#define Cr 1.0 // LambdaT = Cr * T^2
#define xi02 1.0 // correlation length squared

// heat conductivity
PRECISION lambdaT(PRECISION T){
    return Cr * T * T;
}

// heat capacity density
PRECISION Cp(PRECISION s, PRECISION rhob, PRECISION corrL2){
    return (s * s / rhob) * (corrL2 / xi02);
}

// correlation length
PRECISION xi(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// derivative of log(phi0) with respect to energy
PRECISION dlnPhi0de(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// derivative of log(xi) with respect to energy
PRECISION dlnXide(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// derivative of log(phi0) with respect to baryon density
PRECISION dlnPhi0drhob(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

// derivative of log(xi) with respect to baryon density
PRECISION dlnXidrhob(PRECISION e, PRECISION rhob){
    return 1.0; // to be fixed
}

/**************************************************************************************************************************************************/
/* Conductivity, heat capacity, correlation length and relaxation coefficients of different slow modes, and their derivatives
/**************************************************************************************************************************************************/

// universal function
PRECISION f2(PRECISION x){
    return 1.0 / (1.0 + x*x); // Eq. (93)
}

// relaxation coefficents of fluctuations, without (Q*xi)f2(Q*xi), just 2*lambdaT/(Cp*xi^2).
PRECISION relaxationCoefficientPhi(PRECISION rhob, PRECISION s, PRECISION T, PRECISION corrL2)
{
    PRECISION lambdat = lambdaT(T);
    PRECISION cp = Cp(s, rhob, corrL2);
    
    return 2 * lambdat/(cp * corrL2);
}

// relaxation coefficents of fluctuations, only work for f2 defined above.
PRECISION relaxationCoefficientPhiQ(PRECISION gammaPhi, PRECISION corrL2, PRECISION Q)
{
    PRECISION Q2 = Q * Q;
    PRECISION qL2 = corrL2 * Q2;
    PRECISION qL4 = qL2 * qL2;
    
    return gammaPhi * (qL2 + qL4);
}


/**************************************************************************************************************************************************/
/* slow modes out and at equilibrium with different Q
/**************************************************************************************************************************************************/

// slow modes with zero Q
PRECISION equilibriumPhi0(PRECISION rhob, PRECISION s, PRECISION corrL2)
{
    return Cp(s, rhob, corrL2) / (rhob * rhob); // slow modes at equilibrium with Q = 0, Eq.(90)
}

// slow modes with nonzero Q
PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION s, PRECISION Q)
{
    PRECISION corrL = xi(e, rhob);
    PRECISION corrL2 = corrL * corrL;
    PRECISION qL = Q * corrL;
    
    return equilibriumPhi0(rhob, s, corrL2) * f2(qL); // Magnitude of mode Q at Equilibrium, Eq. (89)
}

// initialization of slow modes
void setInitialConditionSlowModes(void * latticeParams, void * hydroParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
#ifdef HydroPlus
    printf("Hydro+ is on...\n");
    
    // initialization of Q vectors
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        Qvec[n] = 0.0 + n * dQvec;
    }

    // initialization of slow mdoes at/out of equilibrium
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                PRECISION es = e[s];
                PRECISION rhobs = rhob[s];
                PRECISION ps = p[s];
                PRECISION Ts = T[s];
                PRECISION alphaBs = muB[s];
                PRECISION entropy = equilibriumEntropy(es, rhobs, ps, Ts, alphaBs);
                
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    
                    PRECISION equiPhiQ = equilibriumPhiQ(es, rhobs, entropy, Qvec[n]);
                    
                    eqPhiQ->phiQ[n][s] = equiPhiQ;
                    eqPhiQp->phiQ[n][s] = equiPhiQ;
                    q->phiQ[n][s] = equiPhiQ;
                    Q->phiQ[n][s] = equiPhiQ;
                }
            }
        }
    }
#endif
}

/**************************************************************************************************************************************************/
/* contributions from slow modes to inferred variables for a specific cell
/**************************************************************************************************************************************************/

// note: the integrands of alpha and beta only work for f2 defined above.
// this function takes e/p/rhob/T/alphaB and slow modes PhiQ/eqPhiQ, then returns variables with contributions from slow modes, including p/T/alphaB
void getPrimaryVariablesFromSlowModes(PRECISION * const __restrict__ p, PRECISION * const __restrict__ T, PRECISION * const __restrict__ alphaB, const PRECISION * const __restrict__ equiPhiQ, const PRECISION * const __restrict__ PhiQ, PRECISION ePrev, PRECISION rhobPrev, PRECISION pPrev, PRECISION TPrev, PRECISION alphaBPrev)
{

    PRECISION corrL = xi(ePrev, rhobPrev); // correlation length
    PRECISION corrL2 = corrL * corrL;
    
    PRECISION s = equilibriumEntropy(ePrev, rhobPrev, pPrev, TPrev, alphaBPrev);
    PRECISION heatC = Cp(s, rhobPrev, corrL2); // heat capacity
    
    PRECISION dlnPhi0_de = dlnPhi0de(ePrev, rhobPrev);
    PRECISION dlnPhi0_drhob = dlnPhi0drhob(ePrev, rhobPrev);
    PRECISION dlnXi_de = dlnXide(ePrev, rhobPrev);
    PRECISION dlnXi_drhob = dlnXidrhob(ePrev, rhobPrev);
    
    PRECISION entropy = 0.0;
    PRECISION alpha = 0.0;
    PRECISION beta = 0.0;
    
    // dQ/(2*pi)^2
    PRECISION facQ = dQvec/(4 * M_PI * M_PI);
    
    // contributions from slow modes to alpha, beta and entropy
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        // ln(Phi/eqPhi) and (Phi/eqPhi-1)
        PRECISION phiRatio = PhiQ[n] / equiPhiQ[n];
        PRECISION phiRatioLog = log(phiRatio);
        PRECISION phiRatioOne = phiRatio - 1;

        // (Q*xi)f2(Q*xi)
        PRECISION qL = Qvec[n] * corrL;
        PRECISION qL2= qL * qL;
        PRECISION qLf2 = qL / (1 + qL2);
        
        // Q^2
        PRECISION Q2 = Qvec[n] * Qvec[n];
        // Q^2*(Phi/eqPhi-1)
        PRECISION QphiRatioOne = Q2 * phiRatioOne;
        // Q^2*ln(Phi/eqPhi)
        PRECISION QphiRatioLog = Q2 * log(phiRatio);
        
        // delta beta, Eq.(106)
        PRECISION intBeta = QphiRatioOne * (dlnPhi0_de - 2 * qLf2 * dlnXi_de);
        beta += intBeta;
        
        // delta alpha, Eq.(106)
        PRECISION intAlpha = QphiRatioOne * (dlnPhi0_drhob - 2 * qLf2 * dlnXi_drhob);
        alpha += intAlpha;
        
        // delta entropy, Eq.(85)
        PRECISION intEntropy = QphiRatioLog - QphiRatioOne;
        entropy += intEntropy;
    }
    
    // contributions from slow modes to entrop, inverse temperature, chemical potential over temperature and pressure
    PRECISION deltaS = facQ * entropy;
    PRECISION deltaAlphaB = - facQ * alpha;
    PRECISION deltaBeta = facQ * beta;
    
    // variables(+) with contribution from slow modes
    *T = 1 / (1/TPrev + deltaBeta);
    *alphaB = alphaBPrev + deltaAlphaB;
    
    PRECISION deltaP = (*T) * (deltaS - (ePrev + pPrev) * deltaBeta + rhobPrev * deltaAlphaB);
    *p = pPrev + deltaP;
}
