#include "fcc.h"

/* -----------------------------------------------------------------------------
 * Constructor, initialize variables for fcc
 * ---------------------------------------------------------------------------*/
FCC::FCC()
{
  NumNei = 12;
  N1stN  = 12;

  memory->create(pos_ideal, NumNei+1, 3, "FCC:pos_ideal");

  pos_ideal[ 0][0] = 0.5;
  pos_ideal[ 0][1] = 0.5;
  pos_ideal[ 0][2] = 0.;
  
  pos_ideal[ 1][0] =-0.5;
  pos_ideal[ 1][1] =-0.5;
  pos_ideal[ 1][2] = 0.;

  pos_ideal[ 2][0] = 0.5;
  pos_ideal[ 2][1] =-0.5;
  pos_ideal[ 2][2] = 0.;
  
  pos_ideal[ 3][0] =-0.5;
  pos_ideal[ 3][1] = 0.5;
  pos_ideal[ 3][2] = 0.;
  
  pos_ideal[ 4][0] = 0.;
  pos_ideal[ 4][1] = 0.5;
  pos_ideal[ 4][2] = 0.5;
  
  pos_ideal[ 5][0] = 0.;
  pos_ideal[ 5][1] =-0.5;
  pos_ideal[ 5][2] =-0.5;

  pos_ideal[ 6][0] =-0.5;
  pos_ideal[ 6][1] = 0.0;
  pos_ideal[ 6][2] = 0.5;
  
  pos_ideal[ 7][0] = 0.5;
  pos_ideal[ 7][1] = 0.0;
  pos_ideal[ 7][2] =-0.5;
  
  pos_ideal[ 8][0] = 0.0;
  pos_ideal[ 8][1] =-0.5;
  pos_ideal[ 8][2] = 0.5;
  
  pos_ideal[ 9][0] = 0.0;
  pos_ideal[ 9][1] = 0.5;
  pos_ideal[ 9][2] =-0.5;
  
  pos_ideal[10][0] = 0.5;
  pos_ideal[10][1] = 0.0;
  pos_ideal[10][2] = 0.5;
  
  pos_ideal[11][0] =-0.5;
  pos_ideal[11][1] = 0.0;
  pos_ideal[11][2] =-0.5;
  
  pos_ideal[12][0] = 0.;
  pos_ideal[12][1] = 0.;
  pos_ideal[12][2] = 0.;
  
return;
}

/* -----------------------------------------------------------------------------
 * method to align local cluster according to fcc
 * ---------------------------------------------------------------------------*/
void FCC::align(double pos_loc[][3])
{
  int nangl = NumNei*(NumNei-1)/2;
  int npair = NumNei/2, counted[NumNei];
  int angends[nangl][2], pairs[npair][2];
  double Ang[nangl], vlen[NumNei];

  // compute the angles between each pair of atoms
  int iangl = 0;
  for (int i=0; i<NumNei; i++){
    vlen[i] = norm2(pos_loc[i]);
    for (int j=i+1; j<NumNei; j++){
      double cos = angle(pos_loc[i],pos_loc[j]);

      Ang[iangl] = cos;
      angends[iangl][0] = i;
      angends[iangl][1] = j;
      iangl++;
    }
  }
  // sort all angles
  for (int i=0; i<iangl; i++)
  for (int j=i+1; j<iangl; j++){
    if (Ang[j] < Ang[i]){
      double atmp = Ang[i];
      Ang[i] = Ang[j];
      Ang[j] = atmp;

      int k1 = angends[i][0];
      int k2 = angends[i][1];
      angends[i][0] = angends[j][0];
      angends[i][1] = angends[j][1];
      angends[j][0] = k1;
      angends[j][1] = k2;
    }
  }

  // pairs with 180 angles are set as pair
  for (int i=0; i<NumNei; i++) counted[i] = 0;
  int ipair = 0;
  for (int i=0; i<iangl; i++){
    int id = angends[i][0], jd = angends[i][1];
    if (counted[id] || counted[jd]) continue;
    pairs[ipair][0] = id; counted[id] = 1;
    pairs[ipair][1] = jd; counted[jd] = 1;
    ipair++;
    if (ipair >= npair) break;
  }

  // store original position
  double pos_init[NumNei][3];
  for (int i=0; i<NumNei; i++)
  for (int j=0; j<3; j++) pos_init[i][j] = pos_loc[i][j];

  // to assign pairs
  int pair_used[npair], index[NumNei];
  for (int i=0; i<npair; i++) pair_used[i] = 0;

  int pmax, pmin;
  int id, jd, iid, jjd, iiid, jjjd;
  // pair most parallel to [110] set as {110}
  double angmax = -1., vij[3], v110[3];
  v110[0] = 1.; v110[1] = 1.; v110[2] = 0.;
  for (int ip=0; ip<npair; ip++){
    id = pairs[ip][0]; jd = pairs[ip][1];
    for (int j=0; j<3; j++) vij[j] = pos_init[jd][j] - pos_init[id][j];
    double ang = angle(vij, v110); ang *= ang;
    if (ang > angmax){
      angmax = ang;
      pmax = ip;
    }
  }
  pair_used[pmax] = 1;
  id = pairs[pmax][0]; jd = pairs[pmax][1];
  for (int j=0; j<3; j++) vij[j] = pos_init[id][j];
  angmax = angle(vij, v110);
  if ( angmax > 0. ){
    index[0] = id; index[1] = jd;
  } else {
    index[0] = jd; index[1] = id;
  }
  id = index[0]; jd = index[1];

  // pair most perpendicular to {110} set as {-110}
  double angmin = 1;
  for (int j=0; j<3; j++) v110[j] = pos_init[id][j] - pos_init[jd][j];
  for (int ip=0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    iid = pairs[ip][0]; jjd = pairs[ip][1];
    for (int j=0; j<3; j++) vij[j] = pos_init[jjd][j] - pos_init[iid][j];
    double ang = angle(vij, v110); ang *= ang;
    if (ang < angmin){
      angmin = ang; pmin = ip;
    }
  }
  pair_used[pmin] = 1;
  iid = pairs[pmin][0]; jjd = pairs[pmin][1];
  if (pos_init[iid][0] > pos_init[jjd][0]){
    index[2] = iid; index[3] = jjd;
  } else {
    index[2] = jjd; index[3] = iid;
  }
  iid = index[2]; jjd = index[3];

  double vm110[3], v001[3]; // define vectors [-110] and [001]
  for (int j=0; j<3; j++) vm110[j] = pos_init[iid][j] - pos_init[jjd][j];
  cross(&v001[0], &vm110[0], &v110[0]);

  // atom with max sum-cosine to 0, 3 and [001] as 011; its counterpart as 0-1-1
  angmax = -1.;
  for (int ip=0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    iiid = pairs[ip][0]; jjjd = pairs[ip][1];
    double ang = angle(pos_init[iiid], pos_init[id]) + angle(pos_init[iiid], pos_init[jjd]) + angle(pos_init[iiid], v001);
    if (ang > angmax){
      angmax = ang;
      pmax = ip; pmin = iiid;
    }
    ang = angle(pos_init[jjjd], pos_init[id]) + angle(pos_init[jjjd], pos_init[jjd]) + angle(pos_init[jjjd], v001);
    if (ang > angmax){
      angmax = ang;
      pmax = ip; pmin = jjjd;
    }
  }
  pair_used[pmax] = 1;
  index[4] = pmin;
  if (pairs[pmax][0] == pmin) index[5] = pairs[pmax][1];
  else index[5] = pairs[pmax][0];

  // atom with max sum-cosine to 1, 3 and [001] as -101; its counterpart as 10-1
  angmax = -1.;
  for (int ip=0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    iiid = pairs[ip][0]; jjjd = pairs[ip][1];
    double ang = angle(pos_init[iiid], pos_init[jd]) + angle(pos_init[iiid], pos_init[jjd]) + angle(pos_init[iiid], v001);
    if (ang > angmax){
      angmax = ang;
      pmax = ip; pmin = iiid;
    }
    ang = angle(pos_init[jjjd], pos_init[jd]) + angle(pos_init[jjjd], pos_init[jjd]) + angle(pos_init[jjjd], v001);
    if (ang > angmax){
      angmax = ang;
      pmax = ip; pmin = jjjd;
    }
  }
  pair_used[pmax] = 1;
  index[6] = pmin;
  if (pairs[pmax][0] == pmin) index[7] = pairs[pmax][1];
  else index[7] = pairs[pmax][0];

  // pair most perpendicular to 011 as {0-11}
  double v011[3];
  for (int j=0; j<3; j++) v011[j] = pos_init[index[4]][j] - pos_init[index[5]][j];
  angmin = 1;
  for (int ip=0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    iiid = pairs[ip][0]; jjjd = pairs[ip][1];
    for (int j=0; j<3; j++) vij[j] = pos_init[jjjd][j] - pos_init[iiid][j];
    double ang = angle(vij, v011); ang *= ang;
    if (ang < angmin){
      angmin = ang; pmin = ip;
    }
  }
  pair_used[pmin] = 1;
  iiid = pairs[pmin][0]; jjjd = pairs[pmin][1];
  angmax = angle(pos_init[iiid], v001);
  if (angmax > 0.){
    index[8] = iiid; index[9] = jjjd;
  } else {
    index[8] = jjjd; index[9] = iiid;
  }

  // the remaining pair as {101}
  for (int ip=0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    pmin = ip;
  }
  pair_used[pmin] = 1;
  iiid = pairs[pmin][0]; jjjd = pairs[pmin][1];
  angmax = angle(pos_init[iiid], v001);
  if (angmax > 0.){
    index[10] = iiid; index[11] = jjjd;
  } else {
    index[10] = jjjd; index[11] = iiid;
  }

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
double FCC::comp_alat(const double vn, const double v0)
{
  return pow(3.*vn+v0, 1./3.);
}

/* -------------------------------------------------------------------------- */
