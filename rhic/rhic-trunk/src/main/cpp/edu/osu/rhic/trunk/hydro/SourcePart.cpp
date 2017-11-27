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
 * Read in the Source terms generated by the particles
/*********************************************************************************************************/
void readInSourcePart(void * latticeParams, void * initCondParams)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;
    //t component
    
    //Initialize Part array
    string line;
    ifstream sourcefile("testsource.txt");
    assert(sourcefile.is_open());
    
//    while (!sourcefile.eof()){
        for(int i = 2; i < nx+2; ++i){
            for(int j = 2; j < ny+2; ++j){
                for(int k = 2; k < nz+2; ++k){
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    if(getline(sourcefile, line,' '))//delimitation character
                    {
                       float fline=atof(line.c_str());
                       Part[s] = (PRECISION) fline;
                       cout<< "Part[ijks]"<< i <<" " << j<<" "  << k <<" " << s <<" " << Part[s] <<endl;
                    }
                }
            }
        }
 //   }
    sourcefile.close();
}

void noSourcePart(void * latticeParams, void * initCondParams)
{
        printf("noSourcePart is started\n");
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
     printf("nx, ny, nz, s=%d, %d, %d, %d\n", nx, ny, nz, (nx+1)+(nx+4)*(ny+1+(ny+4)*(nz+1)));
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                Part[s] = 0;
                cout<< "Part[ijks]"<< i <<" " << j<<" "  << k <<" " << s <<" " << Part[s] <<endl;
            }
        }
        printf("a\n");
    }
    printf("noSourcePart is done\n");
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
	switch (particleSourceType)
    {
		case 0: {
			printf("read in particle source\n");
            readInSourcePart(latticeParams, initCondParams);
			return;
		        }
		case 1: {
			printf("no particle source...\n");
            noSourcePart(latticeParams, initCondParams);
            printf("b\n");
		        }
	}
}
