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

#define dQvec 10.0 // difference between Q vectors of slow modes
#define Q0 0.0

/**************************************************************************************************************************************************/
/* Conductivity, heat capacity, correlation length and relaxation coefficients of different slow modes, and their derivatives
/**************************************************************************************************************************************************/

// for testing
#define Cr 1.0 // LambdaT = Cr * T^2
#define xi0 1.0
#define xi02 1.0 // correlation length squared
#define sigmae 50.0
#define sigman 10.0
#define ec 30.0
#define rhobc 30.0
#define Tc 200.0 // to be fixed

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
    return xi0 * exp(- deltaE * deltaE/(2*sigmae*sigmae) - deltaRhob * deltaRhob/(2*sigman*sigman)) + xi0; // to be fixed
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
    printf("Hydro+ is on, number of slow modes is %d, Q0 is %f, dQ is %f...\n",NUMBER_SLOW_MODES, Q0, dQvec);
    
    // initialization of Q vectors
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        Qvec[n] = Q0 + n * dQvec;
    }

    // initialization of slow mdoes at/out of equilibrium
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                PRECISION es = e[s];
                PRECISION rhobs = rhob[s];
                PRECISION seqs = seq[s];
                
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    
                    PRECISION equiPhiQ = equilibriumPhiQ(es, rhobs, seqs, Qvec[n]);
                    
                    eqPhiQ->phiQ[n][s] = equiPhiQ;
                    q->phiQ[n][s] = equiPhiQ;
                    
                    //printf("Initialization: equiPhiQ = %f\n",equiPhiQ);
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
void getPressurePlusFromSlowModes(PRECISION * const __restrict__ pPlus, const PRECISION * const __restrict__ equiPhiQ, const PRECISION * const __restrict__ PhiQ, PRECISION eIn, PRECISION rhobIn, PRECISION pIn, PRECISION TIn, PRECISION alphaBIn, PRECISION sIn)
{

    PRECISION corrL = xi(eIn, rhobIn); // correlation length
    PRECISION corrL2 = corrL * corrL;
    
    PRECISION heatC = Cp(sIn, rhobIn, corrL2); // heat capacity
    
    // for testing
    PRECISION dlnXi_de = dlnXide(eIn, rhobIn);
    PRECISION dlnXi_drhob = dlnXidrhob(eIn, rhobIn);
    PRECISION dlnPhi0_de = dlnPhi0de(TIn, sIn, dlnXi_de);
    PRECISION dlnPhi0_drhob = dlnPhi0drhob(alphaBIn, rhobIn, sIn, dlnXi_drhob);
    
    PRECISION entropy = 0.0;
    PRECISION alpha = 0.0;
    PRECISION beta = 0.0;
    
    // dQ/(2*pi)^2
    PRECISION facQ = dQvec/(4 * M_PI * M_PI);
    
    // contributions from slow modes to alpha, beta and entropy
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        // ln(Phi/eqPhi) and (Phi/eqPhi-1)
        PRECISION phiRatio = PhiQ[n] / (equiPhiQ[n] + 1e-15);
        PRECISION phiRatioLog = log(phiRatio);
        PRECISION phiRatioOne = phiRatio - 1;
        //printf("PhiQ[n]=%f\t equiPhiQ[n]=%f\t phiRatio=%f\t phiRatioLog=%f\n",PhiQ[n],equiPhiQ[n],phiRatio,phiRatioLog);

        // (Q*xi)f2(Q*xi)
        PRECISION Q = Qvec[n] + 0.5 * dQvec;
        PRECISION qL = Q * corrL;
        PRECISION qL2= qL * qL;
        PRECISION qLf2 = qL / (1 + qL2);
        
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
    }
    
    // contributions from slow modes to entrop, inverse temperature, chemical potential over temperature and pressure
    PRECISION deltaS = facQ * entropy;
    PRECISION deltaAlphaB = - facQ * alpha;
    PRECISION deltaBeta = facQ * beta;
    
    //printf("deltaS=%f,\t deltaAlphaB=%f,\t deltaBeta=%f.\n",deltaS,deltaAlphaB,deltaBeta);
    
    // variables(+) with contribution from slow modes
    PRECISION T = 1 / (1/TIn + deltaBeta);
    
    PRECISION deltaP = T * (deltaS - (eIn + pIn) * deltaBeta + rhobIn * deltaAlphaB);
    *pPlus = pIn;// + deltaP;
    //printf("deltaP=%f\t p+=%f\n",deltaP,pIn+deltaP);
}
