//
//  HydroAnalysis.cpp
//  
//
//  Created by Lipei Du on 10/25/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/PrimaryVariables.h"
#include "../include/EquationOfState.h"
#include "../include/HydroAnalysis.h"
#include "../include/HydroPlus.h"


void outputAnalysis(double t, const char *outputDir, void * latticeParams)
{
    FILE *fp;
    char fname[255];
    sprintf(fname, "%s/AnalysisData.dat", outputDir);
    fp=fopen(fname, "a+");
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double x,y,z;
    
    int i,j,k;
    int s;
    
    /*double v2t,becc,eecc,v2t1,v2t2;
    v2t = 0;
    v2t1 = 0;
    v2t2 = 0;
    becc = 0;
    eecc = 0;
    double bymx = 0;
    double bxy = 0;
    double bypx = 0;
    double eymx = 0;
    double exy = 0;
    double eypx = 0;*/
    
    PRECISION phiQ;
    
    //k=(nz+3)/2;
    //j=(ny+3)/2;
    //i=(nx+3)/2;
    for(k = 2; k < nz+2; ++k) {
        z = (k-2 - (nz-1)/2.)*dz;
        for(j = 2; j < ny+2; ++j) {
            y = (j-2 - (ny-1)/2.)*dy;
            for(i = 2; i < nx+2; ++i) {
                x = (i-2 - (nx-1)/2.)*dx;
                s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                //double tt=Ttt(e[s],p[s],u->ut[s],q->pitt[s]);
                //double tx=Ttx(e[s],p[s],u->ut[s],u->ux[s],q->pitx[s]);
                //double ty=Tty(e[s],p[s],u->ut[s],u->uy[s],q->pity[s]);
                //double tn=Ttn(e[s],p[s],u->ut[s],u->un[s],q->pitn[s]);
/*#ifndef PIMUNU
                double pixx=0;
                double piyy=0;
#else
                double pixx=q->pixx[s];
                double piyy=q->piyy[s];
#endif
                double xx=Txx(e[s],p[s],u->ux[s],pixx);
                //double xy=Txy(e[s],p[s],u->ux[s],u->uy[s],q->pixy[s]);
                //double xn=Txn(e[s],p[s],u->ux[s],u->un[s],q->pixn[s]);
                double yy=Tyy(e[s],p[s],u->uy[s],piyy);
                //double yn=Tyn(e[s],p[s],u->uy[s],u->un[s],q->piyn[s]);
                //double nn=Tnn(e[s],p[s],u->un[s],q->pinn[s],t);
                
                bymx = bymx + (y*y - x*x)*rhob[s];
                bxy = bxy + x*y*rhob[s];
                bypx = bypx + (y*y + x*x)*rhob[s];
                eymx = eymx + (y*y - x*x)*e[s];
                exy = exy + x*y*e[s];
                eypx = eypx + (y*y + x*x)*e[s];
                v2t1 = v2t1 + (xx-yy);
                v2t2 = v2t2 + (xx+yy);*/
                
                //if(x==0&&y==0)
               // phiQ = q->phiQ[0][s];
                
                PRECISION eIn = e[s];
                PRECISION rhobIn = rhob[s];
                PRECISION pIn = p[s];
                PRECISION TIn = T[s];
                PRECISION alphaBIn = alphaB[s];
                PRECISION equiPhiQ[NUMBER_SLOW_MODES];
                PRECISION PhiQ[NUMBER_SLOW_MODES];
                
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    
                    equiPhiQ[n] = eqPhiQ->phiQ[n][s];
                    PhiQ[n] = q->phiQ[n][s];//for test
                }
                
                PRECISION p,T,alphaB;

                getPrimaryVariablesFromSlowModes(&p, &T, &alphaB, equiPhiQ, PhiQ, eIn, rhobIn, pIn, TIn, alphaBIn);
            }
        }
    }
    
    /*becc = sqrt(bymx*bymx + 4*bxy*bxy)/bypx;
    eecc = sqrt(eymx*eymx + 4*exy*exy)/eypx;
    v2t = v2t1/v2t2;
    
    fprintf(fp, "%.3f\t%.8f\t%.8f\t%.8f\n",t,v2t,becc,eecc);*/
    
    //fprintf(fp, "%.3f\t%.8f\n",t,phiQ);
    exit(0);
    
    fclose(fp);
}


