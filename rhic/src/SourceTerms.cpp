/*
 * SourceTerms.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include <stdio.h> // for printf
#include <math.h> // for math functions

#include "../include/SourceTerms.h"
#include "../include/FiniteDifference.h"
#include "../include/EnergyMomentumTensor.h"
#include "../include/DynamicalVariables.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h" // for const params
#include "../include/EquationOfState.h" // for bulk terms
#include "../include/DynamicalSources.h"//Lipei
#include "../include/FluxLimiter.h"//Lipei

//#define USE_CARTESIAN_COORDINATES

//parameters for the analytic parameterization of the bulk viscosity \zeta/S
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

//Transport coefficents of the baryon evolution; Lipei
#define Cb 0.4

inline PRECISION bulkViscosityToEntropyDensity(PRECISION T) {
	PRECISION x = T/1.01355;
	if(x > 1.05)
		return LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
	else if(x < 0.995)
		return LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
	else
		return A_1*x*x + A_2*x - A_3;
}

const PRECISION delta_pipi = 1.33333;
const PRECISION tau_pipi = 1.42857;
const PRECISION delta_PiPi = 0.666667;
const PRECISION lambda_piPi = 1.2;

inline PRECISION baryonDiffusionCoefficient(PRECISION T, PRECISION rhob, PRECISION mub, PRECISION e, PRECISION p){
    PRECISION HyCotangent = 1/tanh(mub);
    if(isnan(HyCotangent))
        printf("kappaB is nan. e=%4e,\t rhob=%4e,\t mub_over_T=%4e,\t T=%4e. \n",e,rhob, mub, T);
    
    return Cb/T * rhob * (0.3333333*HyCotangent - rhob*T/(e+p));
}

void setDissipativeSourceTerms(PRECISION * const __restrict__ pimunuRHS, PRECISION * const __restrict__ nbmuRHS,
        PRECISION nbt, PRECISION nbx, PRECISION nby, PRECISION nbn, PRECISION rhob, PRECISION mub,
        PRECISION Nablat_alphaB, PRECISION Nablax_alphaB, PRECISION Nablay_alphaB, PRECISION Nablan_alphaB,
        PRECISION T, PRECISION t, PRECISION e, PRECISION p,
		PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un, PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp,
		PRECISION pitt, PRECISION pitx, PRECISION pity, PRECISION pitn, PRECISION pixx, PRECISION pixy,
        PRECISION pixn, PRECISION piyy, PRECISION piyn, PRECISION pinn, PRECISION Pi,
		PRECISION dxut, PRECISION dyut, PRECISION dnut, PRECISION dxux, PRECISION dyux, PRECISION dnux,
		PRECISION dxuy, PRECISION dyuy, PRECISION dnuy, PRECISION dxun, PRECISION dyun, PRECISION dnun,
        PRECISION dkvk, PRECISION d_etabar, PRECISION d_dt
) {
	//*********************************************************\
	//* Temperature dependent shear transport coefficients
	//*********************************************************/
	PRECISION taupiInv = T / 5  / d_etabar; //d_etabar=shearViscosityToEntropyDensity
	PRECISION beta_pi = (e + p) / 5;

	//*********************************************************\
	//* Temperature dependent bulk transport coefficients
	//*********************************************************/
	PRECISION cs2 = speedOfSoundSquared(e, rhob);
	PRECISION a = 1.0/3.0 - cs2;
	PRECISION a2 = a*a;
	PRECISION beta_Pi = 15*a2*(e+p);
	PRECISION lambda_Pipi = 8*a/5;

    PRECISION zetabar = bulkViscosityToEntropyDensity(T);
	PRECISION tauPiInv = 15*a2*T/zetabar;

	PRECISION ut2 = ut * ut;
	PRECISION un2 = un * un;
	PRECISION t2 = t * t;
	PRECISION t3 = t * t2;

	// time derivatives of u
	PRECISION dtut = (ut - utp) / d_dt;
	PRECISION dtux = (ux - uxp) / d_dt;
	PRECISION dtuy = (uy - uyp) / d_dt;
	PRECISION dtun = (un - unp) / d_dt;

	//*********************************************************\
	//* covariant derivatives
	//*********************************************************/
	PRECISION Dut = ut*dtut + ux*dxut + uy*dyut + un*dnut + t*un*un;
	PRECISION DuxUpper = ut*dtux + ux*dxux + uy*dyux + un*dnux;
	PRECISION Dux = -DuxUpper;
	PRECISION DuyUpper = ut*dtuy + ux*dxuy + uy*dyuy + un*dnuy;
	PRECISION Duy = -DuyUpper;
	PRECISION DunUpper = ut*dtun + ux*dxun + uy*dyun + un*dnun + 2*ut*un/t;
	PRECISION Dun = -t2*DunUpper;

	PRECISION dut = Dut -t*un*un;
	PRECISION dux = ut*dtux + ux*dxux + uy*dyux + un*dnux;
	PRECISION duy = ut*dtuy + ux*dxuy + uy*dyuy + un*dnuy;
	PRECISION dun = ut*dtun + ux*dxun + uy*dyun + un*dnun;

	//*********************************************************\
	//* expansion rate
	//*********************************************************/
	PRECISION theta = ut / t + dtut + dxux + dyuy + dnun;

	//*********************************************************\
	//* shear tensor
	//*********************************************************/
	PRECISION stt = -t * ut * un2 + (dtut - ut * dut) + (ut2 - 1) * theta / 3;
	PRECISION stx = -(t * un2 * ux) / 2 + (dtux - dxut) / 2 - (ux * dut + ut * dux) / 2 + ut * ux * theta / 3;
	PRECISION sty = -(t * un2 * uy) / 2 + (dtuy - dyut) / 2 - (uy * dut + ut * duy) / 2 + ut * uy * theta / 3;
	PRECISION stn = -un * (2 * ut2 + t2 * un2) / (2 * t) + (dtun - dnut / t2) / 2 - (un * dut + ut * dun) / 2 + ut * un * theta / 3;
	PRECISION sxx = -(dxux + ux * dux) + (1 + ux*ux) * theta / 3;
	PRECISION sxy = -(dxuy + dyux) / 2 - (uy * dux + ux * duy) / 2	+ ux * uy * theta / 3;
	PRECISION sxn = -ut * ux * un / t - (dxun + dnux / t2) / 2 - (un * dux + ux * dun) / 2 + ux * un * theta / 3;
	PRECISION syy = -(dyuy + uy * duy) + (1 + uy*uy) * theta / 3;
	PRECISION syn = -ut * uy * un / t - (dyun + dnuy / t2) / 2 - (un * duy + uy * dun) / 2 + uy * un * theta / 3;
	PRECISION snn = -ut * (1 + 2 * t2 * un2) / t3 - dnun / t2 - un * dun + (1 / t2 + un2) * theta / 3;

	//*********************************************************\
	//* vorticity tensor
	//*********************************************************/
    PRECISION wtx = 0;//(dtux + dxut) / 2 + (ux * dut - ut * dux) / 2 + t * un2 * ux / 2;
	PRECISION wty = 0;//(dtuy + dyut) / 2 + (uy * dut - ut * duy) / 2 + t * un2 * uy / 2;
	PRECISION wtn = 0;//(t2 * dtun + 2 * t * un + dnut) / 2 + (t2 * un * dut - ut * Dun) + t3 * un*un2 / 2;
	PRECISION wxy = 0;//(dyux - dxuy) / 2 + (uy * dux - ux * duy) / 2;
	PRECISION wxn = 0;//(dnux - t2 * dxun) / 2 + (t2 * un * dux - ux * Dun) / 2;
	PRECISION wyn = 0;//(dnuy - t2 * dyun) / 2 + (t2 * un * duy - uy * Dun) / 2;
	// anti-symmetric vorticity components
	PRECISION wxt = 0;//wtx;
	PRECISION wyt = 0;//wty;
	PRECISION wnt = 0;//wtn / t2;
	PRECISION wyx = 0;//-wxy;
	PRECISION wnx = 0;//-wxn / t2;
	PRECISION wny = 0;//-wyn / t2;

#ifdef PIMUNU
	//*********************************************************\
	//* I1
	//*********************************************************/
	PRECISION I1tt = 2 * ut * (pitt * Dut + pitx * Dux + pity * Duy + pitn * Dun);
	PRECISION I1tx = (pitt * ux + pitx * ut) * Dut + (pitx * ux + pixx * ut) * Dux + (pity * ux + pixy * ut) * Duy + (pitn * ux + pixn * ut) * Dun;
	PRECISION I1ty = (pitt * uy + pity * ut) * Dut + (pitx * uy + pixy * ut) * Dux + (pity * uy + piyy * ut) * Duy + (pitn * uy + piyn * ut) * Dun;
	PRECISION I1tn = (pitt * un + pitn * ut) * Dut + (pitx * un + pixn * ut) * Dux + (pity * un + piyn * ut) * Duy + (pitn * un + pinn * ut) * Dun;
	PRECISION I1xx = 2 * ux * (pitx * Dut + pixx * Dux + pixy * Duy + pixn * Dun);
	PRECISION I1xy = (pitx * uy + pity * ux) * Dut + (pixx * uy + pixy * ux) * Dux + (pixy * uy + piyy * ux) * Duy + (pixn * uy + piyn * ux) * Dun;
	PRECISION I1xn = (pitx * un + pitn * ux) * Dut + (pixx * un + pixn * ux) * Dux + (pixy * un + piyn * ux) * Duy + (pixn * un + pinn * ux) * Dun;
	PRECISION I1yy = 2 * uy * (pity * Dut + pixy * Dux + piyy * Duy + piyn * Dun);
	PRECISION I1yn = (pity * un + pitn * uy) * Dut + (pixy * un + pixn * uy) * Dux + (piyy * un + piyn * uy) * Duy + (piyn * un + pinn * uy) * Dun;
	PRECISION I1nn = 2 * un * (pitn * Dut + pixn * Dux + piyn * Duy + pinn * Dun);

	//*********************************************************\
	//* I2
	//*********************************************************/
	PRECISION I2tt = theta * pitt;
	PRECISION I2tx = theta * pitx;
	PRECISION I2ty = theta * pity;
	PRECISION I2tn = theta * pitn;
	PRECISION I2xx = theta * pixx;
	PRECISION I2xy = theta * pixy;
	PRECISION I2xn = theta * pixn;
	PRECISION I2yy = theta * piyy;
	PRECISION I2yn = theta * piyn;
	PRECISION I2nn = theta * pinn;


	//*********************************************************\
	//* I3
	//*********************************************************/
	PRECISION I3tt = 2 * (pitx * wtx + pity * wty + pitn * wtn);
	PRECISION I3tx = pitt * wxt + pity * wxy + pitn * wxn + pixx * wtx + pixy * wty + pixn * wtn;
	PRECISION I3ty = pitt * wyt + pitx * wyx + pitn * wyn + pixy * wtx + piyy * wty + piyn * wtn;
	PRECISION I3tn = pitt * wnt + pitx * wnx + pity * wny + pixn * wtx + piyn * wty + pinn * wtn;
	PRECISION I3xx = 2 * (pitx * wxt + pixy * wxy + pixn * wxn);
	PRECISION I3xy = pitx * wyt + pity * wxt + pixx * wyx + piyy * wxy + pixn * wyn + piyn * wxn;
	PRECISION I3xn = pitx * wnt + pitn * wxt + pixx * wnx + pixy * wny + piyn * wxy + pinn * wxn;
	PRECISION I3yy = 2 * (pity * wyt + pixy * wyx + piyn * wyn);
	PRECISION I3yn = pity * wnt + pitn * wyt + pixy * wnx + pixn * wyx + piyy * wny + pinn * wyn;
	PRECISION I3nn = 2 * (pitn * wnt + pixn * wnx + piyn * wny);

	//*********************************************************\
	//* I4
	//*********************************************************/
	PRECISION ux2 = ux*ux;
	PRECISION uy2 = uy*uy;
	PRECISION ps  = pitt*stt-2*pitx*stx-2*pity*sty+pixx*sxx+2*pixy*sxy+piyy*syy - 2*pitn*stn*t2+2*pixn*sxn*t2+2*piyn*syn*t2+pinn*snn*t2*t2;

	PRECISION I4tt = (pitt * stt - pitx * stx  - pity * sty - t2 * pitn * stn) - (1 - ut2) * ps / 3;
	PRECISION I4tx = (pitt * stx + pitx * stt) / 2 - (pitx * sxx + pixx * stx) / 2 - (pity * sxy + pixy * sty) / 2 - t2 * (pitn * sxn + pixn * stn) / 2 + (ut * ux) * ps / 3;
	PRECISION I4ty = (pitt * sty + pity * stt) / 2 - (pitx * sxy + pixy * stx) / 2 - (pity * syy + piyy * sty) / 2 - t2 * (pitn * syn + piyn * stn) / 2 + (ut * uy) * ps / 3;
	PRECISION I4tn = (pitt * stn + pitn * stt) / 2 - (pitx * sxn + pixn * stx) / 2 - (pity * syn + piyn * sty) / 2 - t2 * (pitn * snn + pinn * stn) / 2 + (ut * un) * ps / 3;
	PRECISION I4xx = (pitx * stx - pixx * sxx  - pixy * sxy - t2 * pixn * sxn) + (1 + ux2) * ps / 3;
	PRECISION I4xy = (pitx * sty + pity * stx) / 2 - (pixx * sxy + pixy * sxx) / 2 - (pixy * syy + piyy * sxy) / 2 - t2 * (pixn * syn + piyn * sxn) / 2 + (ux * uy) * ps / 3;
	PRECISION I4xn = (pitx * stn + pitn * stx) / 2 - (pixx * sxn + pixn * sxx) / 2 - (pixy * syn + piyn * sxy) / 2 - t2 * (pixn * snn + pinn * sxn) / 2 + (ux * un) * ps / 3;
	PRECISION I4yy = (pity * sty - pixy * sxy  - piyy * syy - t2 * piyn * syn) + (1 + uy2) * ps / 3;
	PRECISION I4yn = (pity * stn + pitn * sty) / 2 - (pixy * sxn + pixn * sxy) / 2 - (piyy * syn + piyn * syy) / 2 - t2 * (piyn * snn + pinn * syn) / 2 + (uy * un) * ps / 3;
	PRECISION I4nn = (pitn * stn - pixn * sxn  - piyn * syn - t2 * pinn * snn) + (1 / t2 + un2) * ps / 3;

	//*********************************************************\
	//* I
	//*********************************************************/
	PRECISION Itt = I1tt + delta_pipi * I2tt - I3tt + tau_pipi * I4tt - lambda_piPi * Pi * stt;
	PRECISION Itx = I1tx + delta_pipi * I2tx - I3tx + tau_pipi * I4tx - lambda_piPi * Pi * stx;
	PRECISION Ity = I1ty + delta_pipi * I2ty - I3ty + tau_pipi * I4ty - lambda_piPi * Pi * sty;
	PRECISION Itn = I1tn + delta_pipi * I2tn - I3tn + tau_pipi * I4tn - lambda_piPi * Pi * stn;
	PRECISION Ixx = I1xx + delta_pipi * I2xx - I3xx + tau_pipi * I4xx - lambda_piPi * Pi * sxx;
	PRECISION Ixy = I1xy + delta_pipi * I2xy - I3xy + tau_pipi * I4xy - lambda_piPi * Pi * sxy;
	PRECISION Ixn = I1xn + delta_pipi * I2xn - I3xn + tau_pipi * I4xn - lambda_piPi * Pi * sxn;
	PRECISION Iyy = I1yy + delta_pipi * I2yy - I3yy + tau_pipi * I4yy - lambda_piPi * Pi * syy;
	PRECISION Iyn = I1yn + delta_pipi * I2yn - I3yn + tau_pipi * I4yn - lambda_piPi * Pi * syn;
	PRECISION Inn = I1nn + delta_pipi * I2nn - I3nn + tau_pipi * I4nn - lambda_piPi * Pi * snn;

	//*********************************************************\
	//* shear stress tensor source terms, i.e. terms on RHS
	//*********************************************************/
	PRECISION dpitt = 2 * beta_pi * stt - pitt * taupiInv - Itt - 2 * un * t * pitn;
	PRECISION dpitx = 2 * beta_pi * stx - pitx * taupiInv - Itx - un * t * pixn;
	PRECISION dpity = 2 * beta_pi * sty - pity * taupiInv - Ity - un * t * piyn;
	PRECISION dpitn = 2 * beta_pi * stn - pitn * taupiInv - Itn - un * t * pinn - (ut * pitn + un * pitt) / t;
	PRECISION dpixx = 2 * beta_pi * sxx - pixx * taupiInv - Ixx;
	PRECISION dpixy = 2 * beta_pi * sxy - pixy * taupiInv - Ixy;
	PRECISION dpixn = 2 * beta_pi * sxn - pixn * taupiInv - Ixn - (ut * pixn + un * pitx) / t;
	PRECISION dpiyy = 2 * beta_pi * syy - piyy * taupiInv - Iyy;
	PRECISION dpiyn = 2 * beta_pi * syn - piyn * taupiInv - Iyn - (ut * piyn + un * pity) / t;
	PRECISION dpinn = 2 * beta_pi * snn - pinn * taupiInv - Inn - 2 * (ut * pinn + un * pitn) / t;
#endif
#ifdef PI
	//*********************************************************\
	//* bulk viscous pressure source terms, i.e. terms on RHS
	//*********************************************************/
	PRECISION dPi = -beta_Pi*theta - Pi*tauPiInv - delta_PiPi*Pi*theta + lambda_Pipi*ps;
#endif
	//*********************************************************\
	//* time derivative of the dissipative quantities
	//*********************************************************/
#ifdef PIMUNU
	pimunuRHS[0] = dpitt / ut + pitt * dkvk;
	pimunuRHS[1] = dpitx / ut + pitx * dkvk;
	pimunuRHS[2] = dpity / ut + pity * dkvk;
	pimunuRHS[3] = dpitn / ut + pitn * dkvk;
	pimunuRHS[4] = dpixx / ut + pixx * dkvk;
	pimunuRHS[5] = dpixy / ut + pixy * dkvk;
	pimunuRHS[6] = dpixn / ut + pixn * dkvk;
	pimunuRHS[7] = dpiyy / ut + piyy * dkvk;
	pimunuRHS[8] = dpiyn / ut + piyn * dkvk;
	pimunuRHS[9] = dpinn / ut + pinn * dkvk;
#endif
#ifdef PI
	pimunuRHS[10] = dPi / ut + Pi * dkvk;
#endif

    //*********************************************************\
    //* for the diffusion current of baryon, by Lipei
    //*********************************************************/
    
#ifdef VMU
 
    PRECISION kappaB = baryonDiffusionConstant(T, mub*T)*T;//baryonDiffusionCoefficient(T, rhob, mub, e, p);//
    PRECISION tau_n = Cb/T;
    PRECISION delta_nn = tau_n;
    PRECISION lambda_nn = 0.60 * tau_n;

    PRECISION facNBI1 = nbt * Dut + nbx * Dux + nby * Duy + nbn * Dun;
    PRECISION NBI1t = ut * facNBI1;
    PRECISION NBI1x = ux * facNBI1;
    PRECISION NBI1y = uy * facNBI1;
    PRECISION NBI1n = un * facNBI1;
    
    PRECISION facNBI2 = 1/tau_n * delta_nn * theta;
    PRECISION NBI2t = facNBI2 * nbt;
    PRECISION NBI2x = facNBI2 * nbx;
    PRECISION NBI2y = facNBI2 * nby;
    PRECISION NBI2n = facNBI2 * nbn;
    
    PRECISION NBI3t = -nbx * wxt - nby * wyt - t2 * nbn * wnt;
    PRECISION NBI3x = -nbt * wtx + nby * wyx + t2 * nbn * wnx;
    PRECISION NBI3y = -nbt * wty + nbx * wxy + t2 * nbn * wny;
    PRECISION NBI3n = 1/t2 * (-nbt * wtn + nbx * wxn + nby * wyn);
    
    PRECISION facNBI4 = lambda_nn/tau_n;
    PRECISION NBI4t = facNBI4 * (stt * nbt - stx * nbx - sty * nby - t2 * stn * nbn);
    PRECISION NBI4x = facNBI4 * (stx * nbt - sxx * nbx - sxy * nby - t2 * sxn * nbn);
    PRECISION NBI4y = facNBI4 * (sty * nbt - sxy * nbx - syy * nby - t2 * syn * nbn);
    PRECISION NBI4n = facNBI4 * (stn * nbt - sxn * nbx - syn * nby - t2 * snn * nbn);
    
    PRECISION GBt   = -t * un/ut * nbn;
    PRECISION GBn   = -1/t * (nbn + un/ut * nbt);

    nbmuRHS[0] = -1/ut * (1/tau_n * nbt - 1/tau_n * kappaB * Nablat_alphaB + NBI1t + NBI2t + NBI3t + NBI4t) + nbt * dkvk + GBt;
    nbmuRHS[1] = -1/ut * (1/tau_n * nbx - 1/tau_n * kappaB * Nablax_alphaB + NBI1x + NBI2x + NBI3x + NBI4x) + nbx * dkvk;
    nbmuRHS[2] = -1/ut * (1/tau_n * nby - 1/tau_n * kappaB * Nablay_alphaB + NBI1y + NBI2y + NBI3y + NBI4y) + nby * dkvk;
    nbmuRHS[3] = -1/ut * (1/tau_n * nbn - 1/tau_n * kappaB * Nablan_alphaB + NBI1n + NBI2n + NBI3n + NBI4n) + nbn * dkvk + GBn;
#endif
}


//***************************************************************************************************************************************************/
void loadSourceTermsX(const PRECISION * const __restrict__ I, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dx
) {
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================
	PRECISION facX = 1/d_dx/2;
	int ptr = 20; // 5 * n (with n = 4 corresponding to pitt)
#ifdef PIMUNU
	PRECISION dxpitt = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpitx = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpity = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpitn = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpixx = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpixy = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpixn = (*(I + ptr + 3) - *(I + ptr + 1)) *facX; ptr+=20;
#endif

	PRECISION ut = u->ut[s];
	PRECISION ux = u->ux[s];

	//=========================================================
	// set dx terms in the source terms of T^\tau\mu
	//=========================================================
	PRECISION vx = ux / ut;
#ifdef PIMUNU
#ifndef PI
	S[0] = dxpitt*vx - dxpitx;
	S[1] = dxpitx*vx - dxpixx;
#else
	PRECISION dxPi = (*(I + ptr + 3) - *(I + ptr + 1)) *facX; ptr+=5;
	S[0] = dxpitt*vx - dxpitx - vx*dxPi;
	S[1] = dxpitx*vx - dxpixx - dxPi;
#endif
	S[2] = dxpity*vx - dxpixy;
	S[3] = dxpitn*vx - dxpixn;
#else
    S[0] = 0.0;
    S[1] = 0.0;
    S[2] = 0.0;
    S[3] = 0.0;
#endif
    
    //=========================================================
    // spatial derivatives of the conserved variables \nb^{\mu}; by Lipei
    //=========================================================
#ifdef VMU
    ptr+=5; // increased by 5 to start with nbt
    PRECISION dxnbt = (*(I + ptr + 3) - *(I + ptr + 1)) *facX; ptr+=5;
    PRECISION dxnbx = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;
#endif
    //=========================================================
    // set dx terms in the source terms of baryon current; by Lipei
    //=========================================================
#ifdef NBMU
#ifdef VMU
    S[NUMBER_CONSERVED_VARIABLES] = dxnbt*vx - dxnbx;
#else
    S[NUMBER_CONSERVED_VARIABLES] = 0.0;
#endif
#endif
}

void loadSourceTermsY(const PRECISION * const __restrict__ J, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dy
) {
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================
	PRECISION facY = 1/d_dy/2;
	int ptr = 20; // 5 * n (with n = 4 corresponding to pitt)
#ifdef PIMUNU
	PRECISION dypitt = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypitx = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypity = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypitn = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=10;
	PRECISION dypixy = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=10;
	PRECISION dypiyy = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypiyn = (*(J + ptr + 3) - *(J + ptr + 1)) *facY; ptr+=10;
#endif

	PRECISION ut = u->ut[s];
	PRECISION uy = u->uy[s];

	//=========================================================
	// set dy terms in the source terms
	//=========================================================
	PRECISION vy = uy / ut;
#ifdef PIMUNU
#ifndef PI
	S[0] = dypitt*vy - dypity;
	S[2] = dypity*vy - dypiyy;
#else
	PRECISION dyPi = (*(J + ptr + 3) - *(J + ptr + 1)) *facY; ptr+=5;
	S[0] = dypitt*vy - dypity - vy*dyPi;
	S[2] = dypity*vy - dypiyy - dyPi;
#endif
	S[1] = dypitx*vy - dypixy;
	S[3] = dypitn*vy - dypiyn;
#else
    S[0] = 0.0;
    S[1] = 0.0;
    S[2] = 0.0;
    S[3] = 0.0;
#endif
    
    //=========================================================
    // spatial derivatives of the conserved variables \nb^{\mu}; by Lipei
    //=========================================================
#ifdef VMU
    ptr+=5;
    PRECISION dynbt = (*(J + ptr + 3) - *(J + ptr + 1)) *facY; ptr+=10;
    PRECISION dynby = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;
#endif
    //=========================================================
    // set dy terms in the source terms of baryon current; by Lipei
    //=========================================================
#ifdef NBMU
#ifdef VMU
    S[NUMBER_CONSERVED_VARIABLES] = dynbt*vy - dynby;
#else
    S[NUMBER_CONSERVED_VARIABLES] = 0.0;
#endif
#endif
}

void loadSourceTermsZ(const PRECISION * const __restrict__ K, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION t,
PRECISION d_dz
) {
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================
	PRECISION facZ = 1/d_dz/2;
	int ptr = 20; // 5 * n (with n = 4 corresponding to pitt)
#ifdef PIMUNU
	PRECISION dnpitt = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpitx = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpity = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpitn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=15;
	PRECISION dnpixn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=10;
	PRECISION dnpiyn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpinn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ; ptr+=5;
#endif

	PRECISION ut = u->ut[s];
	PRECISION un = u->un[s];

	//=========================================================
	// set dn terms in the source terms
	//=========================================================
	PRECISION vn = un / ut;
#ifdef PIMUNU
#ifndef PI
	S[0] = dnpitt*vn - dnpitn;
	S[3] = dnpitn*vn - dnpinn;
#else
	PRECISION dnPi = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ; ptr+=5;
	S[0] = dnpitt*vn - dnpitn - vn*dnPi;
	S[3] = dnpitn*vn - dnpinn - dnPi/pow(t,2);
#endif
	S[1] = dnpitx*vn - dnpixn;
	S[2] = dnpity*vn - dnpiyn;
#else
    S[0] = 0.0;
    S[1] = 0.0;
    S[2] = 0.0;
    S[3] = 0.0;
#endif
    
    //=========================================================
    // spatial derivatives of the conserved variables \nb^{\mu}; by Lipei
    //=========================================================
#ifdef VMU
    ptr+=5;
    PRECISION dnnbt = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ; ptr+=15;
    PRECISION dnnbn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;
#endif
    //=========================================================
    // set dn terms in the source terms of baryon current; by Lipei
    //=========================================================
#ifdef NBMU
#ifdef VMU
    S[NUMBER_CONSERVED_VARIABLES] = dnnbt*vn - dnnbn;
#else
    S[NUMBER_CONSERVED_VARIABLES] = 0.0;
#endif
#endif
}

void loadSourceTerms2(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u,
PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp,
PRECISION t, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ pvec,
int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_etabar, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz,
const DYNAMICAL_SOURCE * const __restrict__ Source, const PRECISION * const __restrict__ rhobvec, const PRECISION * const __restrict__ muBvec, const PRECISION * const __restrict__ muBp, const PRECISION * const __restrict__ Tvec, PRECISION Tp
) {
	//=========================================================
	// conserved variables
	//=========================================================
	PRECISION ttt = Q[0];
	PRECISION ttx = Q[1];
	PRECISION tty = Q[2];
	PRECISION ttn = Q[3];
#ifdef PIMUNU
	PRECISION pitt = Q[4];
	PRECISION pitx = Q[5];
	PRECISION pity = Q[6];
	PRECISION pitn = Q[7];
	PRECISION pixx = Q[8];
	PRECISION pixy = Q[9];
	PRECISION pixn = Q[10];
	PRECISION piyy = Q[11];
	PRECISION piyn = Q[12];
	PRECISION pinn = Q[13];
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
	PRECISION Pi = Q[14];
#else
	PRECISION Pi = 0;
#endif
#ifdef NBMU
    PRECISION Nbt = Q[NUMBER_CONSERVED_VARIABLES];
#else
    PRECISION Nbt = 0;
#endif
#ifdef VMU
    PRECISION nbt = Q[NUMBER_CONSERVED_VARIABLES+1];
    PRECISION nbx = Q[NUMBER_CONSERVED_VARIABLES+2];
    PRECISION nby = Q[NUMBER_CONSERVED_VARIABLES+3];
    PRECISION nbn = Q[NUMBER_CONSERVED_VARIABLES+4];
#else
    PRECISION nbt = 0;
    PRECISION nbx = 0;
    PRECISION nby = 0;
    PRECISION nbn = 0;
#endif

	//=========================================================
	// primary variables
	//=========================================================
	PRECISION *utvec = u->ut;
	PRECISION *uxvec = u->ux;
	PRECISION *uyvec = u->uy;
	PRECISION *unvec = u->un;

	PRECISION p  = pvec[s];
	PRECISION ut = utvec[s];
	PRECISION ux = uxvec[s];
	PRECISION uy = uyvec[s];
	PRECISION un = unvec[s];

	//=========================================================
	// spatial derivatives of primary variables
	//=========================================================
	PRECISION facX = 1/d_dx/2;
	PRECISION facY = 1/d_dy/2;
	PRECISION facZ = 1/d_dz/2;
    
	// dx of u^{\mu} components
	PRECISION dxut = (*(utvec + s + 1) - *(utvec + s - 1)) * facX;
	PRECISION dxux = (*(uxvec + s + 1) - *(uxvec + s - 1)) * facX;
	PRECISION dxuy = (*(uyvec + s + 1) - *(uyvec + s - 1)) * facX;
	PRECISION dxun = (*(unvec + s + 1) - *(unvec + s - 1)) * facX;
    
    /*PRECISION utxp1 = utvec[s+1];
    PRECISION utxm1 = utvec[s-1];
    PRECISION uxxp1 = uxvec[s+1];
    PRECISION uxxm1 = uxvec[s-1];
    PRECISION uyxp1 = uyvec[s+1];
    PRECISION uyxm1 = uyvec[s-1];
    PRECISION unxp1 = unvec[s+1];
    PRECISION unxm1 = unvec[s-1];
    
    PRECISION dxut = approximateDerivative(utxm1,ut,utxp1) * facX * 2;
    PRECISION dxux = approximateDerivative(uxxm1,ux,uxxp1) * facX * 2;
    PRECISION dxuy = approximateDerivative(uyxm1,uy,uyxp1) * facX * 2;
    PRECISION dxun = approximateDerivative(unxm1,un,unxp1) * facX * 2;*/
    
	// dy of u^{\mu} components
	PRECISION dyut = (*(utvec + s + d_ncx) - *(utvec + s - d_ncx)) * facY;
	PRECISION dyux = (*(uxvec + s + d_ncx) - *(uxvec + s - d_ncx)) * facY;
	PRECISION dyuy = (*(uyvec + s + d_ncx) - *(uyvec + s - d_ncx)) * facY;
	PRECISION dyun = (*(unvec + s + d_ncx) - *(unvec + s - d_ncx)) * facY;
    
    /*PRECISION utyp1 = utvec[s+d_ncx];
    PRECISION utym1 = utvec[s-d_ncx];
    PRECISION uxyp1 = uxvec[s+d_ncx];
    PRECISION uxym1 = uxvec[s-d_ncx];
    PRECISION uyyp1 = uyvec[s+d_ncx];
    PRECISION uyym1 = uyvec[s-d_ncx];
    PRECISION unyp1 = unvec[s+d_ncx];
    PRECISION unym1 = unvec[s-d_ncx];
    
    PRECISION dyut = approximateDerivative(utym1,ut,utyp1) * facY * 2;
    PRECISION dyux = approximateDerivative(uxym1,ux,uxyp1) * facY * 2;
    PRECISION dyuy = approximateDerivative(uyym1,uy,uyyp1) * facY * 2;
    PRECISION dyun = approximateDerivative(unym1,un,unyp1) * facY * 2;*/
    
	// dn of u^{\mu} components
	int stride = d_ncx * d_ncy;
	PRECISION dnut = (*(utvec + s + stride) - *(utvec + s - stride)) * facZ;
	PRECISION dnux = (*(uxvec + s + stride) - *(uxvec + s - stride)) * facZ;
	PRECISION dnuy = (*(uyvec + s + stride) - *(uyvec + s - stride)) * facZ;
	PRECISION dnun = (*(unvec + s + stride) - *(unvec + s - stride)) * facZ;
    
    /*PRECISION utnp1 = utvec[s+stride];
    PRECISION utnm1 = utvec[s-stride];
    PRECISION uxnp1 = uxvec[s+stride];
    PRECISION uxnm1 = uxvec[s-stride];
    PRECISION uynp1 = uyvec[s+stride];
    PRECISION uynm1 = uyvec[s-stride];
    PRECISION unnp1 = unvec[s+stride];
    PRECISION unnm1 = unvec[s-stride];
    
    PRECISION dnut = approximateDerivative(utnm1,ut,utnp1) * facZ * 2;
    PRECISION dnux = approximateDerivative(uxnm1,ux,uxnp1) * facZ * 2;
    PRECISION dnuy = approximateDerivative(uynm1,uy,uynp1) * facZ * 2;
    PRECISION dnun = approximateDerivative(unnm1,un,unnp1) * facZ * 2;*/
    
	// pressure
    PRECISION dxp = (*(pvec + s + 1) - *(pvec + s - 1)) * facX;
	PRECISION dyp = (*(pvec + s + d_ncx) - *(pvec + s - d_ncx)) * facY;
	PRECISION dnp = (*(pvec + s + stride) - *(pvec + s - stride)) * facZ;
    
    /*PRECISION pc = pvec[s];
    PRECISION pxp1 = pvec[s+1];
    PRECISION pxm1 = pvec[s-1];
    PRECISION pyp1 = pvec[s+d_ncx];
    PRECISION pym1 = pvec[s-d_ncx];
    PRECISION pnp1 = pvec[s+stride];
    PRECISION pnm1 = pvec[s-stride];
    
    PRECISION dxp = approximateDerivative(pxm1,pc,pxp1) * facX * 2;
    PRECISION dyp = approximateDerivative(pym1,pc,pyp1) * facY * 2;
    PRECISION dnp = approximateDerivative(pnm1,pc,pnp1) * facZ * 2;*/

	//=========================================================
	// T^{\mu\nu} source terms
	//=========================================================
	PRECISION tnn = Tnn(evec[s],p+Pi,un,pinn,t);
	PRECISION vx = ux/ut;
	PRECISION vy = uy/ut;
	PRECISION vn = un/ut;
	PRECISION dxvx = (dxux - vx * dxut)/ ut;
	PRECISION dyvy = (dyuy - vy * dyut)/ ut;
	PRECISION dnvn = (dnun - vn * dnut)/ ut;
    PRECISION dkvk = dxvx + dyvy + dnvn;
    
    // added dynamical source terms; lipei
	S[0] = Source->sourcet[s] - (ttt / t + t * tnn) + dkvk*(pitt-p-Pi) - vx*dxp - vy*dyp - vn*dnp;
	S[1] = Source->sourcex[s] - ttx/t -dxp + dkvk*pitx;
	S[2] = Source->sourcey[s] - tty/t -dyp + dkvk*pity;
	S[3] = Source->sourcen[s] - 3*ttn/t -dnp/pow(t,2) + dkvk*pitn;
#ifdef USE_CARTESIAN_COORDINATES
	S[0] = Source->sourcet[s] + dkvk*(pitt-p-Pi) - vx*dxp - vy*dyp - vn*dnp;
	S[1] = Source->sourcex[s] - dxp + dkvk*pitx;
	S[2] = Source->sourcey[s] - dyp + dkvk*pity;
	S[3] = Source->sourcen[s] - dnp + dkvk*pitn;
#endif

    //=========================================================
    // N^{\mu} source terms; by Lipei
    //=========================================================
#ifdef NBMU
    S[NUMBER_CONSERVED_VARIABLES] = Source->sourceb[s] - Nbt/t + dkvk*nbt;
#ifdef USE_CARTESIAN_COORDINATES
    S[NUMBER_CONSERVED_VARIABLES] = Source->sourceb[s] + dkvk*nbt;
#endif
#endif

    //Calculate the gradient of chemical potential to temperature ratio
#ifdef VMU
    PRECISION es = evec[s];
    PRECISION T  = Tvec[s];
    PRECISION rhobs = rhobvec[s];
    PRECISION mubs  = muBvec[s];
    PRECISION mubps  = muBp[s];
    
    PRECISION dtmub = (mubs - mubps ) / d_dt;
    PRECISION dxmub = (*(muBvec + s + 1) - *(muBvec + s - 1)) * facX;
    PRECISION dymub = (*(muBvec + s + d_ncx) - *(muBvec + s - d_ncx)) * facY;
    PRECISION dnmub = (*(muBvec + s + stride) - *(muBvec + s - stride)) * facZ;
    
    //termX[s] = baryonDiffusionCoefficient(T, rhobs, mubs, es, pc);//baryonDiffusionConstant(T, mubs*T)*T;//
    
    /*PRECISION mubxp1 = muBvec[s+1];
    PRECISION mubxp2 = muBvec[s+2];
    PRECISION mubxm1 = muBvec[s-1];
    PRECISION mubxm2 = muBvec[s-2];
    PRECISION mubyp1 = muBvec[s+d_ncx];
    PRECISION mubyp2 = muBvec[s+2*d_ncx];
    PRECISION mubym1 = muBvec[s-d_ncx];
    PRECISION mubym2 = muBvec[s-2*d_ncx];
    PRECISION mubnp1 = muBvec[s+stride];
    PRECISION mubnp2 = muBvec[s+2*stride];
    PRECISION mubnm1 = muBvec[s-stride];
    PRECISION mubnm2 = muBvec[s-2*stride];
    
    PRECISION dtmub, dxmub, dymub, dnmub;
    
    if(rhobs<0.02){
        dtmub = (mubs - mubps ) / d_dt;
        dxmub = approximateDerivative(mubxm1,mubs,mubxp1) * facX * 2;
        dymub = approximateDerivative(mubym1,mubs,mubyp1) * facY * 2;
        dnmub = approximateDerivative(mubnm1,mubs,mubnp1) * facZ * 2;
    }else{
        dtmub = (mubs - mubps ) / d_dt;
        dxmub = (mubxm2 - 8*mubxm2 + mubs + 8*mubxp1 - mubxp2) * facX/6;
        dymub = (mubym2 - 8*mubym2 + mubs + 8*mubyp1 - mubyp2) * facY/6;
        dnmub = (mubnm2 - 8*mubnm2 + mubs + 8*mubnp1 - mubnp2) * facZ/6;
    }*/

    PRECISION ukdk_mub = ut * dtmub + ux * dxmub + uy * dymub + un * dnmub;
    
    PRECISION Nablat_alphaB =  dtmub - ut * ukdk_mub;
    PRECISION Nablax_alphaB = -dxmub - ux * ukdk_mub;
    PRECISION Nablay_alphaB = -dymub - uy * ukdk_mub;
    PRECISION Nablan_alphaB = -1/pow(t,2)*dnmub - un * ukdk_mub;

#else
    PRECISION es = evec[s];
    PRECISION T  = Tvec[s];
    PRECISION rhobs = 0;
    PRECISION mubs  = 0;
    
    PRECISION Nablat_alphaB = 0;
    PRECISION Nablax_alphaB = 0;
    PRECISION Nablay_alphaB = 0;
    PRECISION Nablan_alphaB = 0;
#endif
    
	//=========================================================
	// \pi^{\mu\nu} source terms and nb^\mu source terms
	//=========================================================
#ifndef IDEAL
	PRECISION pimunuRHS[NUMBER_DISSIPATIVE_CURRENTS];
    PRECISION nbmuRHS[NUMBER_PROPAGATED_VMU_COMPONENTS];//Lipei
	setDissipativeSourceTerms(pimunuRHS, nbmuRHS, nbt, nbx, nby, nbn, rhobs, mubs, Nablat_alphaB, Nablax_alphaB, Nablay_alphaB, Nablan_alphaB,
            T, t, es, p, ut, ux, uy, un, utp, uxp, uyp, unp, pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn, Pi,
			dxut, dyut, dnut, dxux, dyux, dnux, dxuy, dyuy, dnuy, dxun, dyun, dnun, dkvk, d_etabar, d_dt);
#ifdef PIMUNU
    for(unsigned int n = 0; n < NUMBER_DISSIPATIVE_CURRENTS; ++n) S[n+4] = pimunuRHS[n];
#endif
#ifdef VMU
    for(unsigned int n = 0; n < NUMBER_PROPAGATED_VMU_COMPONENTS; ++n) S[n+1+NUMBER_CONSERVED_VARIABLES] = nbmuRHS[n];//Source terms for baryon diffusion current
#endif
#endif
}
