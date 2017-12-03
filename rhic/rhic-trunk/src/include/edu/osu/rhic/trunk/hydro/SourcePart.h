/*
 * SourcePart.h
 *
 *  Created on: Nov 25, 2017
 *      Author: Lipei
 */

#ifndef SOURCEPART_H_
#define SOURCEPART_H_

void readInSourcePart(void * latticeParams, void * initCondParams);
void noSourcePart(void * latticeParams, void * initCondParams);
void setSourcePart(void * latticeParams, void * initCondParams, void * hydroParams);

#endif /* SOURCEPART_H_ */
