#include "diamond.h"

/* -----------------------------------------------------------------------------
 * Constructor, initialize variables for Diamond
 * ---------------------------------------------------------------------------*/
Diamond::Diamond()
{
  NumNei = 4;
  N1stN  = 4;

  memory->create(pos_ideal, NumNei+1, 3, "Diamond:pos_ideal");

  pos_ideal[ 0][0] = 0.25;
  pos_ideal[ 0][1] = 0.25;
  pos_ideal[ 0][2] = 0.25;
  
  pos_ideal[ 1][0] =-0.25;
  pos_ideal[ 1][1] =-0.25;
  pos_ideal[ 1][2] = 0.25;

  pos_ideal[ 2][0] = 0.25;
  pos_ideal[ 2][1] =-0.25;
  pos_ideal[ 2][2] =-0.25;
  
  pos_ideal[ 3][0] =-0.25;
  pos_ideal[ 3][1] = 0.25;
  pos_ideal[ 3][2] =-0.25;
  
  pos_ideal[ 4][0] = 0.;
  pos_ideal[ 4][1] = 0.;
  pos_ideal[ 4][2] = 0.;
  
return;
}

/* -----------------------------------------------------------------------------
 * method to align local cluster according to diamond
 * ---------------------------------------------------------------------------*/
void Diamond::align(double pos_loc[][3])
{
  // store original position
  double pos_init[NumNei][3];
  for (int i=0; i<NumNei; i++)
  for (int j=0; j<3; j++) pos_init[i][j] = pos_loc[i][j];

  // select first two as {111}  and {-1-11}
  int index[NumNei];
  int id, jd, iid, jjd;
  id = 0; jd = 1; iid = 2; jjd = 3;

  // most parallel to {111}x{-1-11} as 
  double vnorm[3];
  cross(&vnorm[0], pos_init[id], pos_init[jd]);
  double ang1 = angle(vnorm, pos_init[iid]);
  double ang2 = angle(vnorm, pos_init[jjd]);
  if (ang1 < ang2){iid=3; jjd=2;}

  // index is now correct
  index[0] = id;
  index[1] = jd;
  index[2] = iid;
  index[3] = jjd;

  // assign crystal pos
  for (int i=0; i<NumNei; i++){
    id = index[i];
    pos_loc[i][0] = pos_init[id][0];
    pos_loc[i][1] = pos_init[id][1];
    pos_loc[i][2] = pos_init[id][2];
  }

return;
}

/* -----------------------------------------------------------------------------
 * method to compute the local lattice constant
 * ---------------------------------------------------------------------------*/
double Diamond::comp_alat(const double vn, const double v0)
{
  return pow(6.*vn+v0+v0, 1./3.);
}

/* -------------------------------------------------------------------------- */
