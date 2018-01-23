//#include "ToyJetClass.h"
#include "edu/osu/rhic/trunk/eos/EquationOfState.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

class jetParton{
public:
  //contravariant momentum four vector of parton p^{\mu} in milne coords
  double p[4];
  //contravariant position four vector of parton x^{\mu} in milne coords
  double x[4];
  //proper time derivative of four momentum
  double dp_dtau[4];
  //mass of parton
  double mass;

  //update the four momentum of parton; just an euler step
  void updateMomentum(double dtau);
  //update the four position of the parton - NOT CORRECT FIX IT
  void updatePosition(double dtau);
  //get parton energy loss based on fluid variables
  void energyLoss(int nx, int ny, int nz, double dt, double dx, double dy, double dz, double *ut, double* ux, double *uy, double* un, double *e);
};

//update the four momentum of parton; just an euler step
void jetParton::updateMomentum(double dtau)
{
  for (int i = 0; i < 4; i++)
  {
    p[i] = p[i] + (dtau * dp_dtau[i]);
  }
}
//update the four position of the parton - NOT CORRECT FIX IT
void jetParton::updatePosition(double dtau)
{
  x[0] = x[0] + dtau;
  for (int i = 1; i < 4; i++)
  {
    x[i] = x[i] + (dtau * p[i] / p[0]);
  }
}


void jetParton::energyLoss(int nx, int ny, int nz, double dt, double dx, double dy, double dz, double *ut, double* ux, double *uy, double* un, double *e)
{

  int ncx = nx + 4;
  int ncy = ny + 4;

  //get the partons coordinate index
  double xmin = (-1.0) * ((double)(nx - 1) / 2.0) * dx;
  int ix = (int)round((x[1] - xmin) / dx);
  //printf("ix = %d \n", ix);
  double ymin = (-1.0) * ((double)(ny - 1) / 2.0) * dy;
  int iy = (int)round((x[2] - ymin) / dy);
  //printf("iy = %d \n", iy);
  double zmin = (-1.0) * ((double)(nz - 1) / 2.0) * dz;
  int iz = (int)round((x[3] - zmin) / dz);
  //printf("iz = %d \n", iz);
  int s = columnMajorLinearIndex(ix, iy, iz, ncx, ncy);

  //safegaurd in case parton leaves grid!
  if ((ix >= 0) && (ix < nx) && (iy >= 0) && (iy < ny) && (iz >= 0) && (iz < nz))
  {
    //this prescription for drag/energy loss should work for the bjorken flow, may not make sense in general!
    //define a temperature dependent jet-medium 'relaxation time'
    double t_R = 3.0 / (effectiveTemperature(e[s]) * effectiveTemperature(e[s]) * effectiveTemperature(e[s]));

    dp_dtau[0] = (-1.0) * ((p[0] / p[0]) - ut[s]) / t_R;
    dp_dtau[1] = (-1.0) * ((p[1] / p[0]) - ux[s]) / t_R;
    dp_dtau[2] = (-1.0) * ((p[2] / p[0]) - uy[s]) / t_R;
    dp_dtau[3] = (-1.0) * ((p[3] / p[0]) - un[s]) / t_R;

    updateMomentum(dt);
  }

  else
  {
    printf("parton has escaped medium with final position {%f, %f, %f, %f} and momentum {%f, %f, %f, %f} \n", x[0],x[1],x[2],x[3],p[0],p[1],p[2],p[3]);
  }

}
