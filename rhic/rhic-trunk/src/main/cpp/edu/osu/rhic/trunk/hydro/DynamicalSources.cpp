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

/*********************************************************************************************************\
 * Dynamical source terms from the jet traversing the medium
/*********************************************************************************************************/

void setDynamicalSources(void * latticeParams, void * initCondParams, double *dp_dtau, double *pos) //dp_dtau is the jet energy loss, pos is the jet position
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;
	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	PRECISION dummy = 0.0;

    double xmin = (-1.0) * ((double)(nx - 1) / 2.0) * dx;
	double ymin = (-1.0) * ((double)(ny - 1) / 2.0) * dy;
	double zmin = (-1.0) * ((double)(nz - 1) / 2.0) * dz;


	//construct an array of the gaussian smeared jet position
	double smearedPosition[ncx * ncy * ncz];
	double width = 0.1; //width of gaussian smearing

	for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                double x = (double)i * dx + xmin;
                double y = (double)j * dy + ymin;
                double z = (double)k * dz + zmin;
                smearedPosition[s] = exp((-1.0)*(pos[1] - x) * (pos[1] - x) / width) * exp((-1.0)*(pos[2] - y) * (pos[2] - y) / width) * exp((-1.0)*(pos[3] - z) * (pos[3] - z) / width);
            }//k
        }//j
    }//i

    //now multiply the smeared position by energy loss corresponding to vector components
	for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                Source->sourcet[s] = -dp_dtau[0] * smearedPosition[s];
                Source->sourcex[s] = -dp_dtau[1] * smearedPosition[s];
                Source->sourcey[s] = -dp_dtau[2] * smearedPosition[s];
                Source->sourcen[s] = -dp_dtau[3] * smearedPosition[s];
                Source->sourceb[s] = dummy;
            }//k
        }//j
    }//i
}

