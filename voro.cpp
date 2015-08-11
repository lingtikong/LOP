/* ------------------------------------------------------------------
 * class to read the output from voro++; prepare for local order
 * parameter evaluation
 *
 * The voro custom output is expected as
 *  %i %q %v %s %F %f %n
 * ------------------------------------------------------------------ */

#include "voro.h"

/* ------------------------------------------------------------------
 * Constructor as the main driver to read the voro output file
 * ------------------------------------------------------------------ */
Voro::Voro(const char *infile, const double box[3], const int nflag, Crystal *crys) : Neigh(crys)
{
  FILE *fp = fopen(infile, "r");
  if (fp == NULL) return;

  pbc[0] = pbc[1] = pbc[2] = 0;
  for (int i=0; i<3; i++){
    if (box[i] > 0.) pbc[i] = 1;
    prd[i] = box[i]; hprd[i] = 0.5*box[i];
  }
  flag_comp_neigh = nflag;

  int idmax = 0, nguess = 5000;
  natom = 0;

  memory->create(x, nguess,3,"Voro:x");
  memory->create(neig_list, nguess,NumNei,"Voro:neig_list");
  vol.clear();
  vol.resize(nguess,0.);

  while (1){
    fgets(oneline, MAXLINE, fp);
    if (feof(fp)) break;

    int id = atoi(strtok(oneline," \t\n\r\f")); // get atom id
    idmax = MAX(idmax, id);

    if (id >= nguess){
      while (id >= nguess) nguess += 1000;
      memory->grow(x,nguess,3,"Voro:grow:x");
      memory->grow(neig_list,nguess,NumNei,"Voro:grow:neig_list");
      vol.resize(nguess,0.);
    }
    x[id][0] = atof(strtok(NULL," \t\n\r\f")); // get atomic position
    x[id][1] = atof(strtok(NULL," \t\n\r\f"));
    x[id][2] = atof(strtok(NULL," \t\n\r\f"));
    vol[id]  = atof(strtok(NULL," \t\n\r\f")); // volume

    int nvoro = atoi(strtok(NULL," \t\n\r\f"));

    if (flag_comp_neigh == 2 || (flag_comp_neigh==1 && nvoro < NumNei)){
      atom_list.push_back(id);
    } else {
      strtok(NULL," \t\n\r\f"); // skip total voro surface area
      local_prop.clear(); local_list.clear();

      for (int i=0; i<nvoro; i++) local_prop.push_back(atof(strtok(NULL," \t\n\r\f")));
      for (int i=0; i<nvoro; i++) local_list.push_back(atoi(strtok(NULL," \t\n\r\f")));

      sort_list(id);
    }

    natom++;
  }
  fclose(fp);

  flag = idmax - natom;

  if (flag == 0){
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
  }

return;
}

/* ------------------------------------------------------------------
 * To supply the "local lattice constant" of atom id
 * ------------------------------------------------------------------ */
double Voro::get_alat(int id)
{
  return a0[id];
}

/* ------------------------------------------------------------------
 * To free memory
 * ------------------------------------------------------------------ */
Voro::~Voro()
{
  atom_list.clear();
  vol.clear(); a0.clear();
}

/* ------------------------------------------------------------------
 * Private method to compute the neighbor list for atoms with voro
 * neighbors less than NumNei, by searching all atoms
 * ------------------------------------------------------------------ */
void Voro::compute_neigh()
{
  // set cutoff for neighbors as half the minimal box length
  double cutoff = MAX(prd[0], MAX(prd[1], prd[2])) * 0.5;
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
        if (pbc[idim]){
          while (Rij[idim] >  hprd[idim]) Rij[idim] -= prd[idim];
          while (Rij[idim] < -hprd[idim]) Rij[idim] += prd[idim];
        }
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
