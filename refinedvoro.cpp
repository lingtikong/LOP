/* ------------------------------------------------------------------
 * class to read the output from voro calculations of dumpana;
 * prepare for local order parameter evaluation
 * ------------------------------------------------------------------ */

#include "refinedvoro.h"

/* ------------------------------------------------------------------
 * Constructor as the main driver to read the dumpana output file
 * ------------------------------------------------------------------ */
RefinedVoro::RefinedVoro(const char *infile, double *box, const int nflag, Crystal *crys) : Neigh(crys)
{
  FILE *fp = fopen(infile, "r");
  if (fp == NULL) return;
  flag_comp_neigh = nflag;

  fgets(oneline,MAXLINE,fp);
  char *ptr = strtok(oneline, " \n\t\r\f");
  strtok(NULL," \n\t\r\f");
  xlo = atof(strtok(NULL," \n\t\r\f"));
  xhi = atof(strtok(NULL," \n\t\r\f"));
  ylo = atof(strtok(NULL," \n\t\r\f"));
  yhi = atof(strtok(NULL," \n\t\r\f"));
  zlo = atof(strtok(NULL," \n\t\r\f"));
  zhi = atof(strtok(NULL," \n\t\r\f"));
  natom = atoi(strtok(NULL," \n\t\r\f"));
  if (natom < 1) return;

  fbox[0] = box[0] = lx  = xhi-xlo;
  fbox[1] = box[1] = ly  = yhi-ylo;
  fbox[2] = box[2] = lz  = zhi-zlo;
  for (int i=0; i<3; i++) hbox[i] = 0.5*fbox[i];

  memory->create(x, natom+1,3,"RefinedVoro:x");
  memory->create(attyp, natom+1,"RefinedVoro:attyp");
  memory->create(neig_list,natom+1,NumNei,"RefinedVoro:neig_list");
  vol.clear();
  vol.resize(natom+1,0.);

  while (!feof(fp)){
    fgets(oneline, MAXLINE, fp);
    if (count_words(oneline) < 9) continue;

    char *ptr = strtok(oneline," \t\n\r\f");

    int id = atoi(ptr);                      // get atom id
    int ip = atoi(strtok(NULL," \t\n\r\f")); // get atom type

    attyp[id] = ip;

    x[id][0] = atof(strtok(NULL," \t\n\r\f")); // get atomic position
    x[id][1] = atof(strtok(NULL," \t\n\r\f"));
    x[id][2] = atof(strtok(NULL," \t\n\r\f"));
    vol[id]  = atof(strtok(NULL," \t\n\r\f")); // volume

    for (int i=0; i<3; i++) ptr = strtok(NULL," \t\n\r\f");
    int nvoro = atoi(ptr);

    if (flag_comp_neigh == 2 || (flag_comp_neigh==1 && nvoro < NumNei)){
      atom_list.push_back(id);
    } else {
      local_prop.clear(); local_list.clear();

      for (int i=0; i<nvoro; i++) local_list.push_back(atoi(strtok(NULL," \t\n\r\f")));
      for (int i=0; i<nvoro; i++) local_prop.push_back(atof(strtok(NULL," \t\n\r\f")));

      sort_list(id);
    }
  }
  fclose(fp);

  if (atom_list.size() > 0) compute_neigh();

  a0.clear();
  a0.push_back(1.);
  for (int id = 1; id<=natom; id++){
    double vloc = 0.;
    for (int i=0; i<N1stN; i++){
      int jd = neig_list[id][i];
      vloc += vol[jd];
    }
    vloc /= double(N1stN);

    a0.push_back( crystal->comp_alat(vloc, vol[id]) );
  }

  flag = 0;
return;
}

/* ------------------------------------------------------------------
 * To supply the "local lattice constant" of atom id
 * ------------------------------------------------------------------ */
double RefinedVoro::get_alat(int id)
{
  return a0[id];
}

/* ------------------------------------------------------------------
 * To
 * ------------------------------------------------------------------ */
RefinedVoro::~RefinedVoro()
{
  atom_list.clear();
  vol.clear(); a0.clear();
}

/* ------------------------------------------------------------------
 * Private method to compute the neighbor list for atoms with voro
 * neighbors less than NumNei, by searching all atoms
 * ------------------------------------------------------------------ */
void RefinedVoro::compute_neigh()
{
  // set cutoff for neighbors as half the minimal box length
  double cutoff = MIN(lx, MIN(ly, lz)) * 0.499;
  double cutoff2 = cutoff * cutoff;
  cutoff2 = MAX(cutoff2, 50.);

  // to get the neighbor list
  double Rij[3];
  for (int i=0; i<atom_list.size(); i++){
    int id = atom_list[i];
    local_prop.clear(); local_list.clear();

    for (int jd=1; jd<=natom; jd++){
      if (jd == id) continue;

      double rij2 = 0.;
      for (int idim=0; idim<3; idim++){
        Rij[idim] = x[jd][idim] - x[id][idim];
        // apply pbc
        while (Rij[idim] >= hbox[idim]) Rij[idim] -= fbox[idim];
        while (Rij[idim] < -hbox[idim]) Rij[idim] += fbox[idim];
        rij2 += Rij[idim]*Rij[idim];
      }
      if (rij2 >= cutoff2) continue;

      local_prop.push_back(-rij2);
      local_list.push_back(jd);
    }
    int nlocal = local_list.size();
    if (nlocal < NumNei) printf("\nWarning: # of neighbors (%d) computed for %d is < NumNei!\n", nlocal, id);

    sort_list(id);
  }

return;
}
