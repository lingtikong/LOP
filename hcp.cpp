#include "hcp.h"
#include <set>

/* -----------------------------------------------------------------------------
 * Constructor, initialize variables for hcp
 * ---------------------------------------------------------------------------*/
HCP::HCP(const double ca0)
{
  ca = ca0;

  NumNei = 12;
  N1stN  = 12;

  memory->create(pos_ideal, NumNei+1, 3, "HCP:pos_ideal");
  const double hrt3 = 0.5*sqrt(3.);
  const double trt3 = sqrt(3.)/3.;
  const double hca  = ca*0.5;

  pos_ideal[ 0][0] = hrt3;
  pos_ideal[ 0][1] =-0.5;
  pos_ideal[ 0][2] = 0.;
  
  pos_ideal[ 1][0] = hrt3;
  pos_ideal[ 1][1] = 0.5;
  pos_ideal[ 1][2] = 0.;

  pos_ideal[ 2][0] = 0.;
  pos_ideal[ 2][1] = 1.;
  pos_ideal[ 2][2] = 0.;
  
  pos_ideal[ 3][0] =-hrt3;
  pos_ideal[ 3][1] = 0.5;
  pos_ideal[ 3][2] = 0.;
  
  pos_ideal[ 4][0] =-hrt3;
  pos_ideal[ 4][1] =-0.5;
  pos_ideal[ 4][2] = 0.;
  
  pos_ideal[ 5][0] = 0.;
  pos_ideal[ 5][1] =-1.;
  pos_ideal[ 5][2] = 0.;

  pos_ideal[ 6][0] = trt3;
  pos_ideal[ 6][1] = 0.;
  pos_ideal[ 6][2] = hca;
  
  pos_ideal[ 7][0] =-trt3*0.5;
  pos_ideal[ 7][1] = 0.5;
  pos_ideal[ 7][2] = hca;
  
  pos_ideal[ 8][0] =-trt3*0.5;
  pos_ideal[ 8][1] =-0.5;
  pos_ideal[ 8][2] = hca;
  
  pos_ideal[ 9][0] = trt3;
  pos_ideal[ 9][1] = 0.;
  pos_ideal[ 9][2] =-hca;
  
  pos_ideal[10][0] =-trt3*0.5;
  pos_ideal[10][1] = 0.5;
  pos_ideal[10][2] =-hca;
  
  pos_ideal[11][0] =-trt3*0.5;
  pos_ideal[11][1] =-0.5;
  pos_ideal[11][2] =-hca;
  
  pos_ideal[12][0] = 0.;
  pos_ideal[12][1] = 0.;
  pos_ideal[12][2] = 0.;
  
return;
}

/* -----------------------------------------------------------------------------
 * method to align local cluster according to hcp
 * ---------------------------------------------------------------------------*/
void HCP::align(double pos_loc[][3])
{
  int hNumNei = NumNei/2;
  double a2 = 0.;
  for (int i=0; i<NumNei; i++) a2 = MAX(a2, norm2(pos_loc[i]));
  a2 *= 1.2;

  int nnei[NumNei], list[NumNei][NumNei];
  int bonded[NumNei][NumNei];
  for (int i=0; i<NumNei; i++)
  for (int j=i; j<NumNei; j++) bonded[i][j] = 0;

  for (int i=0; i<NumNei; i++) nnei[i] = 0;
  for (int i=0;   i<NumNei; i++)
  for (int j=i+1; j<NumNei; j++){
    double vij[3];
    for (int idim=0; idim<3; idim++) vij[idim] = pos_loc[j][idim]-pos_loc[i][idim];
    double r2 = norm2(vij);
    if (r2 < a2){
      list[i][nnei[i]] = j;
      list[j][nnei[j]] = i;

      bonded[i][j] = bonded[j][i] = 1;
      nnei[i]++; nnei[j]++;
    }
  }

  std::set<int> bondednbs;
  int env[NumNei];
  for (int i=0; i<NumNei; i++){
    bondednbs.clear();
    for (int jj=0;    jj<nnei[i]; jj++)
    for (int kk=jj+1; kk<nnei[i]; kk++){
      int j=list[i][jj], k=list[i][kk];
      if (bonded[j][k]){bondednbs.insert(j); bondednbs.insert(k);}
    }
    env[i] = bondednbs.size();
  }
  bondednbs.clear();

  int ids[NumNei];
  for (int i=0; i<NumNei; i++) ids[i] = i;
  for (int i=0;   i<NumNei; i++)
  for (int j=i+1; j<NumNei; j++){
    if (env[j] < env[i]){
      double swap = env[i];
      env[i] = env[j];
      env[j] = swap;

      int idum = ids[i];
      ids[i] = ids[j];
      ids[j] = idum;
    }
  }
//  for (int i=0; i<NumNei; i++) printf("%d ", env[i]); printf("\n");
  
  double angs[NumNei][NumNei];
  for (int i=0;   i<NumNei; i++)
  for (int j=i+1; j<NumNei; j++) angs[i][j] = angs[j][i] = angle(pos_loc[i], pos_loc[j]);

  int nmax = hNumNei*(hNumNei-1)/2;
  double Ang[nmax];
  int angends[nmax][2];
  int iang = 0;
  for (int ii=0;    ii<hNumNei; ii++)
  for (int jj=ii+1; jj<hNumNei; jj++){
    int i=ids[ii], j=ids[jj];
    Ang[iang] = angs[i][j];
    angends[iang][0] = ii;
    angends[iang][1] = jj;
    iang++;
  }
  for (int ii=0;    ii<iang; ii++)
  for (int jj=ii+1; jj<iang; jj++){
    if (Ang[jj] < Ang[ii]){
      double dtmp = Ang[ii];
      Ang[ii] = Ang[jj];
      Ang[jj] = dtmp;
      int i0 = angends[ii][0];
      int i1 = angends[ii][1];
      angends[ii][0] = angends[jj][0];
      angends[ii][1] = angends[jj][1];
      angends[jj][0] = i0;
      angends[jj][1] = i1;
    }
  }

  // store original position
  double pos_init[NumNei][3];
  for (int i=0; i<NumNei; i++)
  for (int idim=0; idim<3; idim++) pos_init[i][idim] = pos_loc[i][idim];

  // to assign pairs
  int index[NumNei];
  // pairs with 180 angles as 0-3, 1-4, and 2-5
  int counted[NumNei];
  for (int i=0; i<NumNei; i++) counted[i] = 0;

  int i03 = -1, i14 = -1, i25 = -1;
  for (int ia=0; ia<iang; ia++){
    int ii = angends[ia][0], jj = angends[ia][1];
    if (counted[ii] || counted[jj]) continue;
    if (i03 < 0) i03 = ia;
    else if (i14 < 0) i14 = ia;
    else i25 = ia;
    counted[ii] = counted[jj] = 1;
    if (i25 >= 0) break;
  }
  index[0] = ids[angends[i03][0]]; index[3] = ids[angends[i03][1]];

  // for 1-4 pair, set the one with acute angle to 0 as 1
  int i1 = ids[angends[i14][0]], i4 = ids[angends[i14][1]];
  if (angs[i1][index[0]] > angs[i4][index[0]]){
    index[1] = i1; index[4] = i4;
  } else {
    index[1] = i4; index[4] = i1;
  }

  // for 2-5 pair, set the one with acute angle to 0 as 5
  int i2 = ids[angends[i25][0]], i5 = ids[angends[i25][1]];
  if (angs[i2][index[0]] < angs[i5][index[0]]){
    index[2] = i2; index[5] = i5;
  } else {
    index[2] = i5; index[5] = i2;
  }

  // check if rotation needed
  std::set<int> hascomm;
  hascomm.clear();
  for (int ii=0; ii<nnei[index[0]]; ii++) hascomm.insert(list[index[0]][ii]);
  for (int jj=0; jj<nnei[index[1]]; jj++) hascomm.insert(list[index[1]][jj]);
  if (nnei[index[0]]+nnei[index[1]] - (int)hascomm.size() < 2){
    int id5 = index[0];
    for (int i=0; i<5; i++) index[i] = index[i+1];
    index[5] = id5;
  }
      
  // set vector norm to 0-3, 1-4 as basal norm
  double v03[3], v14[3], vnorm[3];
  for (int idim=0; idim<3; idim++){
    v03[idim] = pos_init[index[0]][idim] - pos_init[index[3]][idim];
    v14[idim] = pos_init[index[1]][idim] - pos_init[index[4]][idim];
  }
  cross(&vnorm[0], &v03[0], &v14[0]);

  // compute angle to norm
  double A2N[NumNei];
  for (int ii=hNumNei; ii<NumNei; ii++){
    int i = ids[ii];
    A2N[i] = angle(pos_init[i], vnorm);
  }

  double cond, majo, cosn, cos0, cos1, cos2, cos3, cos4, cos5;
  // atom with smallest sum angle to 0, 1, and norm as 6
  int pmax6;
  double amax6 = -10.;
  for (int ii=hNumNei; ii<NumNei; ii++){
    int i = ids[ii];
    cosn = A2N[i];
    cos0 = angs[i][index[0]];
    cos1 = angs[i][index[1]];
    cos2 = angs[i][index[2]];
    cos3 = angs[i][index[3]];
    cos4 = angs[i][index[4]];
    cos5 = angs[i][index[5]];

    majo = cos0 + cos1 + 3.*cosn;
    cond = 3*majo - cos2*cos2 - cos5*cos5 - cos3 - cos4;
    if (cond > amax6){ amax6 = cond; pmax6 = ii; }
  }
  index[6] = ids[pmax6]; counted[pmax6] = 1;

  // atoms with smallest sum angle to 0, 1 and -norm as 9
  int pmax9;
  double amax9 = -10.;
  for (int ii=hNumNei; ii<NumNei; ii++){
    if (counted[ii]) continue;
    int i = ids[ii];
    cosn = A2N[i];
    cos0 = angs[i][index[0]];
    cos1 = angs[i][index[1]];
    cos2 = angs[i][index[2]];
    cos3 = angs[i][index[3]];
    cos4 = angs[i][index[4]];
    cos5 = angs[i][index[5]];

    majo = cos0 + cos1 - 3.*cosn;
    cond = 3*majo - cos3 - cos4 - cos2*cos2 - cos5*cos5;
    if (cond > amax9){ amax9 = cond; pmax9 = ii; }
  }
  index[9] = ids[pmax9]; counted[pmax9] = 1;

  // atom with smallest sum angle to 2, 3, and norm as 7
  int pmax7;
  double amax7 = -10.;
  for (int ii=hNumNei; ii<NumNei; ii++){
    if (counted[ii]) continue;
    int i=ids[ii];
    cosn = A2N[i];
    cos0 = angs[i][index[0]];
    cos1 = angs[i][index[1]];
    cos2 = angs[i][index[2]];
    cos3 = angs[i][index[3]];
    cos4 = angs[i][index[4]];
    cos5 = angs[i][index[5]];

    majo = cos2 + cos3 + (cosn+cosn);
    cond = 3*majo - cos0 - cos5 - cos1*cos1 - cos4*cos4;
    if (cond > amax7){ amax7 = cond; pmax7 = ii; }
  }
  index[7] = ids[pmax7]; counted[pmax7] = 1;

  // atom with smallest sum angle to 4, 5, and norm as 8
  int pmax8;
  double amax8 = -10.;
  for (int ii=hNumNei; ii<NumNei; ii++){
    if (counted[ii]) continue;
    int i = ids[ii];
    cosn = A2N[i];
    cos0 = angs[i][index[0]];
    cos1 = angs[i][index[1]];
    cos2 = angs[i][index[2]];
    cos3 = angs[i][index[3]];
    cos4 = angs[i][index[4]];
    cos5 = angs[i][index[5]];

    majo = cos4 + cos5 + (cosn+cosn);
    cond = 3*majo - cos1 - cos2 - cos0*cos0 - cos3*cos3;
    if (cond > amax8){ amax8 = cond; pmax8 = ii; }
  }
  index[8] = ids[pmax8]; counted[pmax8] = 1;

  // atom with smallest sum angle to 2, 3, and -norm as 10
  int pmax10;
  double amax10 = -10.;
  for (int ii=hNumNei; ii<NumNei; ii++){
    if (counted[ii]) continue;
    int i=ids[ii];
    cosn = A2N[i];
    cos0 = angs[i][index[0]];
    cos1 = angs[i][index[1]];
    cos2 = angs[i][index[2]];
    cos3 = angs[i][index[3]];
    cos4 = angs[i][index[4]];
    cos5 = angs[i][index[5]];

    majo = cos2 + cos3 - (cosn+cosn);
    cond = 3*majo - cos0 - cos5 - cos1*cos1 - cos4*cos4;
    if (cond > amax10){ amax10 = cond; pmax10 = ii; }
  }
  index[10] = ids[pmax10]; counted[pmax10] = 1;

  // the remaining one as 11
  int pmax11;
  for (int ii=hNumNei; ii<NumNei; ii++){
    if (counted[ii]) continue;
    pmax11 = ii;
  }
  index[11] = ids[pmax11];

  // assign crystal pos
  for (int i=0; i<NumNei; i++){
    int id = index[i];
    pos_loc[i][0] = pos_init[id][0];
    pos_loc[i][1] = pos_init[id][1];
    pos_loc[i][2] = pos_init[id][2];
  }

return;
}

/* -----------------------------------------------------------------------------
 * method to compute the local lattice constant
 * ---------------------------------------------------------------------------*/
double HCP::comp_alat(const double vn, const double v0)
{
  const double fac = 2./(sqrt(3.)*ca);

  return pow((vn+v0)*fac, 1./3.);
}

/* -------------------------------------------------------------------------- */
