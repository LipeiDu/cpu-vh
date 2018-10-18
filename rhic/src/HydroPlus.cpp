//
//  HydroPlus.cpp
//  
//
//  Created by Lipei Du on 10/17/18.
//

// Equation indices from PRD 98 (2018) 036006

#include <math.h>
#include <cmath>

#include "../include/DynamicalVariables.h"
#include "../include/HydroPlus.h"

#define Cp 1.0
#define lambdaT 1.0

SLOW_MODES_Q *Qvec;
SLOW_MODES *phiQ;

//correlation length

PRECISION xi(PRECISION e, PRECISION rhob){
    return 1.0; // correlation length, to be fixed
}

//fluctuations at equilibrium

PRECISION f2(PRECISION x){
    return 1.0 / (1.0 + x*x); // Eq. (93)
}

PRECISION equilibriumPhi0(PRECISION e, PRECISION rhob){
    return Cp / (rhob * rhob); // Cp needs to be fixed, Eq.(90)
}

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
    
    return lambdaT/(Cp * corrL * corrL) * (qL2 + qL4);
}

//source terms of extra modes, Eq. (76)

void setSlowModesSourceTerms(PRECISION * const __restrict__ PhiQRHS, const PRECISION * const __restrict__ equilibriumPhiQ, const PRECISION * const __restrict__ PhiQ, const SLOW_MODES_Q * const __restrict__ Qvec, PRECISION e, PRECISION rhob, PRECISION ut, PRECISION dkvk){
    
    PRECISION gammaQ[NUMBER_SLOW_MODES];
    PRECISION Q[NUMBER_SLOW_MODES];
    PRECISION utInv = 1.0/ut;
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        *Q[n] = Qvec->qvec[n];
        gammaQ[n] = relaxationCoefficientPhiQ(e, rhob, Q[n]);
        PhiQRHS[n] = utInv * (-2) * gammaQ[n] * (PhiQ[n] - equilibriumPhiQ[n]) + PhiQ[n] * dkvk;
    }
    
}

void loadSlowModesSourceTerms(const SLOW_MODES * const __restrict__ equilibriumSlowModes, const SLOW_MODES * const __restrict__ SlowModes, PRECISION * const __restrict__ S, const SLOW_MODES_Q * const __restrict__ Qvec, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ rhobvec, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz){ // Modes, slow modes; S, source terms;
    
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
        
        *equilibriumPhiQvec[n] = equilibriumSlowModes->phiQ[n];
        equilibriumPhiQ[n] = equilibriumPhiQvec[n][s];
        
        *PhiQvec[n] = SlowModes->phiQ[n];
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
