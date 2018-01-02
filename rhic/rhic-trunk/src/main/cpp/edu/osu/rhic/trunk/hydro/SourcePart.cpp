/*
 * SourcePart.c
 *
 *  Created on: Nov 25, 2017
 *  Author: Lipei
 */

#include <math.h> // for math functions
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP
#include <iostream>
#include <istream>
#include <fstream>
#include <cassert>
#include <string>

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/trunk/hydro/SourcePart.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"

using namespace std;
/*********************************************************************************************************\
 * Initialize the Source terms from particles
/*********************************************************************************************************/

void readInSourcePart(void * latticeParams, void * initCondParams)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
    
    FILE *sourcefile;
    
    sourcefile = fopen ("testsource.dat","r");
    if(sourcefile==NULL){
        printf("The source file testsource.dat was not opened...\n");
        exit(-1);
    }
    else
    {
      fseek(sourcefile,0L,SEEK_SET);
      for(int i = 2; i < nx+2; ++i){
         for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
               int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
               fscanf(sourcefile,"%lf %lf %lf %lf %lf", & Part->partt[s], & Part->partx[s], & Part->party[s], & Part->partn[s], & Part->partb[s]);
               printf("%lf %lf %lf %lf %lf\n", Part->partt[s], Part->partx[s], Part->party[s], Part->partn[s], Part->partb[s]);
             }
          }
       }
    }

    fclose(sourcefile);
}

void noSourcePart(void * latticeParams, void * initCondParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                Part->partt[s] = 0;
                Part->partx[s] = 0;
                Part->party[s] = 0;
                Part->partn[s] = 0;
                Part->partb[s] = 0;
            }
        }
    }
}

/*********************************************************************************************************\
 * Source terms from particles
 * 	0 - add these source terms
 *	1 - not add these source terms
/*********************************************************************************************************/

void setSourcePart(void * latticeParams, void * initCondParams, void * hydroParams)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int particleSourceType = initCond->particleSourceType;
	printf("Source terms from particles: ");
	switch (particleSourceType) {
		case 0:{
			printf("read in particle sources\n");
            readInSourcePart(latticeParams, initCondParams);
            printf("Particle sources have been read in.\n");
			return;}
		case 1:{
			printf("set particle sources as 0\n");
            noSourcePart(latticeParams, initCondParams);
            printf("Particle sources initialized.\n");
            return;}
	}
}
