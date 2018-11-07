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

#define dQvec 0.4 // difference between Q vectors of slow modes

/**************************************************************************************************************************************************/
/* Conductivity, heat capacity, correlation length and relaxation coefficients of different slow modes, and their derivatives
/**************************************************************************************************************************************************/

// for testing
#define Cr 1.0 // LambdaT = Cr * T^2
#define xi0 1.0
#define xi02 1.0 // correlation length squared
#define sigmae 50.0
#define sigman 10.0
#define ec 250.0
#define rhobc 40.0

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
    PRECISION deltaE = e - ec;
    PRECISION deltaRhob = rhob - rhobc;
    return xi0 * exp(- deltaE * deltaE/(2*sigmae*sigmae) - deltaRhob * deltaRhob/(2*sigman*sigman)); // to be fixed
}

// derivative of log(xi) with respect to energy
PRECISION dlnXide(PRECISION e, PRECISION rhob){
    return -((e - ec)/(sigmae*sigmae)); // to be fixed
}

// derivative of log(xi) with respect to baryon density
PRECISION dlnXidrhob(PRECISION e, PRECISION rhob){
    return -((rhob - rhobc)/(sigman*sigman)); // to be fixed
}

// derivative of log(phi0) with respect to energy
PRECISION dlnPhi0de(PRECISION T, PRECISION s, PRECISION dlnXi_de){
    return 2/(s*T) + 2 * dlnXi_de; // to be fixed
}

// derivative of log(phi0) with respect to baryon density
PRECISION dlnPhi0drhob(PRECISION alphaB, PRECISION rhob, PRECISION s, PRECISION dlnXi_drhob){
    return -2 * alphaB / s - 3 / rhob + 2 * dlnXi_drhob; // to be fixed
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
                    
                    //printf("equiPhiQ=%f",equiPhiQ);
                    
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
void getPrimaryVariablesFromSlowModes(PRECISION * const __restrict__ p, PRECISION * const __restrict__ T, PRECISION * const __restrict__ alphaB, const PRECISION * const __restrict__ equiPhiQ, const PRECISION * const __restrict__ PhiQ, PRECISION eIn, PRECISION rhobIn, PRECISION pIn, PRECISION TIn, PRECISION alphaBIn)
{

    PRECISION corrL = xi(eIn, rhobIn); // correlation length
    PRECISION corrL2 = corrL * corrL;
    
    //printf("eIn =%f\t rhobIn =%f\n",eIn,rhobIn);
    
    PRECISION s = equilibriumEntropy(eIn, rhobIn, pIn, TIn, alphaBIn);
    PRECISION heatC = Cp(s, rhobIn, corrL2); // heat capacity
    
    // for testing
    PRECISION dlnXi_de = dlnXide(eIn, rhobIn);
    PRECISION dlnXi_drhob = dlnXidrhob(eIn, rhobIn);
    PRECISION dlnPhi0_de = dlnPhi0de(TIn, s, dlnXi_de);
    PRECISION dlnPhi0_drhob = dlnPhi0drhob(alphaBIn, rhobIn, s, dlnXi_drhob);
    
    
    //printf("corrL =%f\t dlnPhi0_de =%f\t dlnPhi0_drhob =%f\t dlnXi_de =%f\t dlnXi_drhob =%f\t \n",corrL,dlnPhi0_de,dlnPhi0_drhob,dlnXi_de,dlnXi_drhob);
    
    PRECISION entropy = 0.0;
    PRECISION alpha = 0.0;
    PRECISION beta = 0.0;
    
    // dQ/(2*pi)^2
    PRECISION facQ = dQvec/(4 * M_PI * M_PI);
    
    // contributions from slow modes to alpha, beta and entropy
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        //printf("n=%d,equiPhiQ[n]=%f\n",n,equiPhiQ[n]);
        
        // ln(Phi/eqPhi) and (Phi/eqPhi-1)
        PRECISION phiRatio = PhiQ[n] / (equiPhiQ[n] + 1e-15);
        PRECISION phiRatioLog = log(phiRatio);
        PRECISION phiRatioOne = phiRatio - 1;

        // (Q*xi)f2(Q*xi)
        PRECISION Q = Qvec[n] + 0.5 * dQvec;
        PRECISION qL = Q * corrL;
        PRECISION qL2= qL * qL;
        PRECISION qLf2 = qL / (1 + qL2);
        
        //printf("phiRatio =%f\t qLf2=%.4e\n",phiRatio,qLf2);
        
        // Q^2
        PRECISION Q2 = Q * Q;
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
        
        //printf("%f\t%f\t%f\t%f\n",Qvec[n],intEntropy,intBeta,intAlpha);
    }
    
    // contributions from slow modes to entrop, inverse temperature, chemical potential over temperature and pressure
    PRECISION deltaS = facQ * entropy;
    PRECISION deltaAlphaB = - facQ * alpha;
    PRECISION deltaBeta = facQ * beta;
    
    //printf("%f\t%f\t%f\n",deltaS,deltaBeta,deltaAlphaB);
    
    // variables(+) with contribution from slow modes
    *T = 1 / (1/TIn + deltaBeta);
    *alphaB = alphaBIn + deltaAlphaB;
    
    PRECISION deltaP = (*T) * (deltaS - (eIn + pIn) * deltaBeta + rhobIn * deltaAlphaB);
    *p = pIn + deltaP;
}
