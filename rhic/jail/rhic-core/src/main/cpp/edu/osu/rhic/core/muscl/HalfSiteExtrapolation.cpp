/*
 * HalfSiteExtrapolation.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/core/muscl/HalfSiteExtrapolation.h"
#include "edu/osu/rhic/core/muscl/FluxLimiter.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

PRECISION rightHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp) {
	return qp - approximateDerivative(q, qp, qpp)/2;
}

PRECISION rightHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp) {
	return q - approximateDerivative(qm, q, qp)/2;
}

PRECISION leftHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp) {
	return q + approximateDerivative(qm, q, qp)/2;
}

PRECISION leftHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp) {
	return qm + approximateDerivative(qmm, qm, q)/2;
}
