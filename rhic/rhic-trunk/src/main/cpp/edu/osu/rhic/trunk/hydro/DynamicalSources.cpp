/*
 * DynamicalSources.cpp
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

#include <iomanip>//by Lipei

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalSources.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"


using namespace std;
//*********************************************************************************************************\
//* Initialize the dynamical source terms
//*********************************************************************************************************/

void readInSource(int n, void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

    // For each time step, the total number of cells in 3D space
    int nElements = nx * ny * nz;
    
    double t0 = hydro->initialProperTimePoint;
    double dt = lattice->latticeSpacingProperTime;
    
    double t = t0 + (n-1)* dt;
    double time;
    
    //printf("t=%lf\n",t);

    FILE *sourcefile;
    char fname[255];
    sprintf(fname, "%s/%s%d.dat", rootDirectory, "../urqmd-source/part2s/output/Sources",n);
    sourcefile = fopen(fname, "r");
    
    FILE *fp;
    char fpname[255];
    sprintf(fpname, "%s/output/sourcet_%d.dat", rootDirectory, n);
    fp=fopen(fpname, "w");
    
    FILE *fpx;
    char fpnamex[255];
    sprintf(fpnamex, "%s/output/sourceb_%d.dat", rootDirectory, n);
    fpx=fopen(fpnamex, "w");

    if(sourcefile==NULL){
        printf("The source file could not be opened...\n");
        exit(-1);
    }
    else
    {
      fseek(sourcefile,0L,SEEK_SET);
        
      //for(int i=0; i<(nElements+1)*(n-1); i++) fscanf(sourcefile,"%*[^\n]%*c");//Skip the title line and all the cells read in by previous steps, (nElements+1) lines

      fscanf(sourcefile,"%*s%le%*c", &time);
      //printf("time=%lf\n",time);
      if(time-t>1.e-20) printf("The dynamical source at a wrong time step is being read in. tSource=%lf, tCode=%lf\n", time, t);
      //if(time==t) printf("The dynamical source starts to be read in at %lf.\n", time);
        
      for(int i = 2; i < nx+2; ++i){
         for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
               int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
               fscanf(sourcefile,"%le %le %le %le %le", & Source->sourcet[s], & Source->sourcex[s], & Source->sourcey[s], & Source->sourcen[s], & Source->sourceb[s]);
               //printf("%le\t %le\t %le\t %le\t %le\n", Source->sourcet[s], Source->sourcex[s], Source->sourcey[s], Source->sourcen[s], Source->sourceb[s]);
               fprintf(fp, "%.3f\t%.3f\t%.3f\t%.8f\n",(i-1 - (201-1)/2) * 0.1,(j-1 - (201-1)/2) * 0.1,(k-1 - (1-1)/2) * 0.1,Source->sourcet[s]);
               fprintf(fpx, "%.3f\t%.3f\t%.3f\t%.8f\n",(i-1 - (201-1)/2) * 0.1,(j-1 - (201-1)/2) * 0.1,(k-1 - (1-1)/2) * 0.1,Source->sourceb[s]);
             }
          }
       }
    }
    
    fclose(fp);//Lipei
    fclose(fpx);
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
                Source->sourcet[s] = 0.0;
                Source->sourcex[s] = 0.0;
                Source->sourcey[s] = 0.0;
                Source->sourcen[s] = 0.0;
                Source->sourceb[s] = 0.0;
            }//k
        }//j
    }//i
}

//*********************************************************************************************************\
//* Dynamical source terms for T^\mu\nu and J^\mu
//* 0 - add these source terms
//*	1 - not add these source terms
/*********************************************************************************************************/

void setSource(int n, void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int sourceType = initCond->sourceType;
    //if ((n-1) % FREQ == 0) printf("(Elapsed time: %.3f ms)\n",delta_time);
	//printf("Dynamical source terms: ");
	switch (sourceType) {
		case 0:{
			//printf("read in dynamical sources\n");
            readInSource(n, latticeParams, initCondParams, hydroParams, rootDirectory);
            //printf("dynamical sources have been read in.\n");
			return;
        }
		case 1:{
			//printf("set dynamical sources as 0\n");
            //noSource(latticeParams, initCondParams);
            //printf("dynamical sources initialized.\n");
            return;
        }
	}
}

//*********************************************************************************************************\
//* Dynamical source terms from the jet traversing the medium
//*********************************************************************************************************/

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

