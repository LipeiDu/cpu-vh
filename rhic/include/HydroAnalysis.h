//
//  HydroAnalysis.h
//  
//
//  Created by Lipei Du on 10/25/18.
//

#ifndef HydroAnalysis_h
#define HydroAnalysis_h

void outputAnalysis(double t, FILE *outputDir, void * latticeParams);
void testEOS();
void testHydroPlus();
void outputGammaQ(double t, const char *pathToOutDir, void * latticeParams);

#endif /* HydroAnalysis_h */
