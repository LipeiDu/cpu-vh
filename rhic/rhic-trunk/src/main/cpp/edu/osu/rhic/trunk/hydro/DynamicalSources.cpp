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
#include "edu/osu/rhic/trunk/hydro/DynamicalSources.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"

using namespace std;
/*********************************************************************************************************\
 * Initialize the dynamical source terms
/*********************************************************************************************************/

void readInSource(void * latticeParams, void * initCondParams)
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
               fscanf(sourcefile,"%lf %lf %lf %lf %lf", & Source->sourcet[s], & Source->sourcex[s], & Source->sourcey[s], & Source->sourcen[s], & Source->sourceb[s]);
               //printf("%lf %lf %lf %lf %lf\n", Source->sourcet[s], Source->sourcex[s], Source->sourcey[s], Source->sourcen[s], Source->sourceb[s]);
             }
          }
       }
    }

    fclose(sourcefile);
}

void noSource(void * latticeParams, void * initCondParams)
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
                Source->sourcet[s] = 0;
                Source->sourcex[s] = 0;
                Source->sourcey[s] = 0;
                Source->sourcen[s] = 0;
                Source->sourceb[s] = 0;
            }//k
        }//j
    }//i
}

/*********************************************************************************************************\
 * Dynamical source terms for T^\mu\nu and J^\mu
 * 	0 - add these source terms
 *	1 - not add these source terms
/*********************************************************************************************************/

void setSource(void * latticeParams, void * initCondParams, void * hydroParams)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int sourceType = initCond->sourceType;
	printf("Dynamical source terms: ");
	switch (sourceType) {
		case 0:{
			printf("read in dynamical sources\n");
            readInSource(latticeParams, initCondParams);
            printf("dynamical sources have been read in.\n");
			return;}
		case 1:{
			printf("set dynamical sources as 0\n");
            noSource(latticeParams, initCondParams);
            printf("dynamical sources initialized.\n");
            return;}
	}
}

void setDynamicalSources(void * latticeParams, void * initCondParams)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	PRECISION dummy = updateDynamicalSources();

	for(int i = 2; i < nx+2; ++i){
			for(int j = 2; j < ny+2; ++j){
					for(int k = 2; k < nz+2; ++k){
							int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
							Source->sourcet[s] = dummy;
							Source->sourcex[s] = dummy;
							Source->sourcey[s] = dummy;
							Source->sourcen[s] = dummy;
							Source->sourceb[s] = dummy;
					}//k
			}//j
	}//i
}

PRECISION updateDynamicalSources()
{
	return 0.0;
}
