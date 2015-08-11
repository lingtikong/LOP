#include "bcc.h"

/* -----------------------------------------------------------------------------
 * Constructor, initialize variables for bcc
 * ---------------------------------------------------------------------------*/
BCC::BCC()
{
  NumNei = 14;
  N1stN  = 8;

  memory->create(pos_ideal, NumNei+1, 3, "BCC:pos_ideal");

  pos_ideal[ 0][0] = 0.;
  pos_ideal[ 0][1] = 0.;
  pos_ideal[ 0][2] = 1.;
  
  pos_ideal[ 1][0] = 0.;
  pos_ideal[ 1][1] = 0.;
  pos_ideal[ 1][2] =-1.;

  pos_ideal[ 2][0] = 1.;
  pos_ideal[ 2][1] = 0.;
  pos_ideal[ 2][2] = 0.;
  
  pos_ideal[ 3][0] =-1.;
  pos_ideal[ 3][1] = 0.;
  pos_ideal[ 3][2] = 0.;
  
  pos_ideal[ 4][0] = 0.;
  pos_ideal[ 4][1] = 1.;
  pos_ideal[ 4][2] = 0.;
  
  pos_ideal[ 5][0] = 0.;
  pos_ideal[ 5][1] =-1.;
  pos_ideal[ 5][2] = 0.;

  pos_ideal[ 6][0] = 0.5;
  pos_ideal[ 6][1] = 0.5;
  pos_ideal[ 6][2] = 0.5;
  
  pos_ideal[ 7][0] =-0.5;
  pos_ideal[ 7][1] =-0.5;
  pos_ideal[ 7][2] =-0.5;
  
  pos_ideal[ 8][0] =-0.5;
  pos_ideal[ 8][1] = 0.5;
  pos_ideal[ 8][2] = 0.5;
  
  pos_ideal[ 9][0] = 0.5;
  pos_ideal[ 9][1] =-0.5;
  pos_ideal[ 9][2] =-0.5;
  
  pos_ideal[10][0] = 0.5;
  pos_ideal[10][1] =-0.5;
  pos_ideal[10][2] = 0.5;
  
  pos_ideal[11][0] =-0.5;
  pos_ideal[11][1] = 0.5;
  pos_ideal[11][2] =-0.5;
  
  pos_ideal[12][0] =-0.5;
  pos_ideal[12][1] =-0.5;
  pos_ideal[12][2] = 0.5;
  
  pos_ideal[13][0] = 0.5;
  pos_ideal[13][1] = 0.5;
  pos_ideal[13][2] =-0.5;
  
  pos_ideal[14][0] = 0.;
  pos_ideal[14][1] = 0.;
  pos_ideal[14][2] = 0.;

return;
}

/* -----------------------------------------------------------------------------
 * method to align local cluster according to bcc
 * ---------------------------------------------------------------------------*/
void BCC::align(double pos_loc[][3])
{
  int nangl = NumNei*(NumNei-1)/2;
  int npair = NumNei/2, counted[NumNei];
  int angends[nangl][2], pairs[npair][2];
  double Ang[nangl], angs[NumNei][NumNei], vlen[NumNei];

  int iangl = 0;
  for (int i=0; i<NumNei; i++){
    vlen[i] = norm2(pos_loc[i]);
    for (int j=i+1; j<NumNei; j++){
      double cos = angle(pos_loc[i],pos_loc[j]);
      angs[i][j] = angs[j][i] = cos;

      Ang[iangl] = cos;
      angends[iangl][0] = i;
      angends[iangl][1] = j;
      iangl++;
    }
  }
  for (int i=0; i<iangl; i++)
  for (int j=i+1; j<iangl; j++){
    if (Ang[j] < Ang[i]){
      double atmp = Ang[i];
      Ang[i] = Ang[j];
      Ang[j] = atmp;

      int k1 = angends[i][0], k2= angends[i][1];
      angends[i][0] = angends[j][0]; angends[i][1] = angends[j][1];
      angends[j][0] = k1; angends[j][1] = k2;
    }
  }

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
  //if (ipair != npair) printf("\nnpair=%d != ipair=%d!\n", npair, ipair);

  // store original position
  double pos_init[NumNei][3];
  for (int i=0; i<NumNei; i++)
  for (int j=0; j<3; j++) pos_init[i][j] = pos_loc[i][j];

  // to assign pairs
  int pair_used[npair], index[NumNei];
  for (int i=0; i<npair; i++) pair_used[i] = 0;

  int pmax, pmin;
  int id, jd, iid, jjd, iiid, jjjd;
  // pair with largest distance + angle as 001
  double lenmax = -1.;
  for (int ip=0; ip<npair; ip++){
    id = pairs[ip][0]; jd = pairs[ip][1];
    double dip = vlen[id] + vlen[jd] - angs[id][jd];
    if (dip > lenmax){
      lenmax = dip;
      pmax = ip;
    }
  }
  pair_used[pmax] = 1;
  id = pairs[pmax][0]; jd = pairs[pmax][1];
  if (pos_init[id][2] > pos_init[jd][2]){
    index[0] = id; index[1] = jd;
  } else {
    index[0] = jd; index[1] = id;
  }
  id = index[0]; jd = index[1];

  // pair with second largest distance as 100
  lenmax = -1.;
  for (int ip=0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    iid = pairs[ip][0]; jjd = pairs[ip][1];
    double dip = vlen[iid] + vlen[jjd] - angs[iid][jjd]
               - (fabs(angs[iid][id]) + fabs(angs[iid][jd]) + fabs(angs[jjd][id]) + fabs(angs[jjd][jd]))*0.5;
    if (dip > lenmax){
      lenmax = dip; pmax = ip;
    }
  }
  pair_used[pmax] = 1;
  iid = pairs[pmax][0]; jjd = pairs[pmax][1];
  if (pos_init[iid][0] > pos_init[jjd][0]){
    index[2] = iid; index[3] = jjd;
  } else {
    index[2] = jjd; index[3] = iid;
  }
  iid = index[2]; jjd = index[3];

  // pair with third largest distance as 010
  lenmax = -1.;
  for (int ip=0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    iiid = pairs[ip][0]; jjjd = pairs[ip][1];
    double dip = vlen[iiid] + vlen[jjjd] - angs[iiid][jjjd] 
               -( fabs(angs[iiid][id]) + fabs(angs[iiid][jd]) + fabs(angs[iiid][iid]) + fabs(angs[iiid][jjd])
               +  fabs(angs[jjjd][id]) + fabs(angs[jjjd][jd]) + fabs(angs[jjjd][iid]) + fabs(angs[jjjd][jjd])) * 0.25;
    if (dip > lenmax){
      lenmax = dip; pmax = ip;
    }
  }
  pair_used[pmax] = 1;
  iiid = pairs[pmax][0]; jjjd = pairs[pmax][1];

  double n010[3],n001[3], n100[3];
  for (int i=0; i<3; i++){ n001[i] = pos_init[id][i]-pos_init[jd][i]; n100[i] = pos_init[iid][i] - pos_init[jjd][i];}
  cross(n010, n001, n100);
  double a0 = angle(n010, pos_init[iiid]), a1 = angle(n010, pos_init[jjjd]);
  if (a0 > a1){
    index[4] = iiid; index[5] = jjjd;
  } else {
    index[4] = jjjd; index[5] = iiid;
  }
  iiid = index[4]; jjjd = index[5];

  //   id: +z;   jd: -z
  //  iid: +x;  jjd: -x
  // iiid: +y; jjjd: -y

  // assign 111
  int ik, jk, flag;
  double cosmax = -3.;
  for (int ip=0; ip < npair; ip++){
    if (pair_used[ip]) continue;
    ik = pairs[ip][0]; jk = pairs[ip][1];
    double cossum = angs[ik][id] + angs[ik][iid] + angs[ik][iiid] + angs[jk][jd] + angs[jk][jjd] + angs[jk][jjjd];
    if (cossum > cosmax){
      flag = 0;
      cosmax = cossum;
      pmax = ip;
    }
    cossum = angs[jk][id] + angs[jk][iid] + angs[jk][iiid] + angs[ik][jd] + angs[ik][jjd] + angs[ik][jjjd];
    if (cossum > cosmax){
      flag = 1;
      cosmax = cossum;
      pmax = ip;
    }
  }
  pair_used[pmax] = 1;
  if (flag == 0){
    index[6] = pairs[pmax][0];
    index[7] = pairs[pmax][1];
  } else {
    index[6] = pairs[pmax][1];
    index[7] = pairs[pmax][0];
  }

  // assign -1-11
  cosmax = -3.;
  for (int ip=0; ip < npair; ip++){
    if (pair_used[ip]) continue;
    ik = pairs[ip][0];
    jk = pairs[ip][1];
    double cossum = angs[ik][jd] + angs[ik][iid] + angs[ik][iiid] + angs[jk][id] + angs[jk][jjd] + angs[jk][jjjd];
    if (cossum > cosmax){
      flag = 0;
      cosmax = cossum;
      pmax = ip;
    }
    cossum = angs[jk][jd] + angs[jk][iid] + angs[jk][iiid] + angs[ik][id] + angs[ik][jjd] + angs[ik][jjjd];
    if (cossum > cosmax){
      flag = 1;
      cosmax = cossum;
      pmax = ip;
    }
  }
  pair_used[pmax] = 1;
  if (flag == 1){
    index[12] = pairs[pmax][0];
    index[13] = pairs[pmax][1];
  } else {
    index[12] = pairs[pmax][1];
    index[13] = pairs[pmax][0];
  }

  // assign -111
  cosmax = -3.;
  for (int ip=0; ip < npair; ip++){
    if (pair_used[ip]) continue;
    ik = pairs[ip][0];
    jk = pairs[ip][1];
    double cossum = angs[ik][id] + angs[ik][jjd] + angs[ik][iiid] + angs[jk][jd] + angs[jk][iid] + angs[jk][jjjd];
    if (cossum > cosmax){
      flag = 0;
      cosmax = cossum;
      pmax = ip;
    }
    cossum = angs[jk][id] + angs[jk][jjd] + angs[jk][iiid] + angs[ik][jd] + angs[ik][iid] + angs[ik][jjjd];
    if (cossum > cosmax){
      flag = 1;
      cosmax = cossum;
      pmax = ip;
    }
  }
  pair_used[pmax] = 1;
  if (flag == 0){
    index[8] = pairs[pmax][0];
    index[9] = pairs[pmax][1];
  } else {
    index[8] = pairs[pmax][1];
    index[9] = pairs[pmax][0];
  }

  // assign 1-11
  for (int ip = 0; ip<npair; ip++){
    if (pair_used[ip]) continue;
    ik = pairs[ip][0];
    jk = pairs[ip][1];
    if (angs[ik][id] > angs[jk][id]){
      index[10] = ik;
      index[11] = jk;
    } else {
      index[10] = jk;
      index[11] = ik;
    }
    break;
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
double BCC::comp_alat(const double vn, const double v0)
{
  return pow(vn+v0, 1./3.);
}

/* -------------------------------------------------------------------------- */
