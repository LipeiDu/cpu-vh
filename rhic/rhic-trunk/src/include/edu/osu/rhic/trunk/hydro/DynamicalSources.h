/*
 * SourcePart.h
 *
 *  Created on: Nov 25, 2017
 *      Author: Lipei
 */

#ifndef SOURCE_H_
#define SOURCE_H_

void readInSource(void * latticeParams, void * initCondParams);
void noSource(void * latticeParams, void * initCondParams);
void setSource(void * latticeParams, void * initCondParams, void * hydroParams);

#endif /* SOURCEPART_H_ */
