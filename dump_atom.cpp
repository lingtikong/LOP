/* ------------------------------------------------------------------
 * class to read dump_atom from LAMMPS; prepare for local order
 * parameter evaluation
 * ------------------------------------------------------------------ */

#include "dump_atom.h"

/* ------------------------------------------------------------------
 * Constructor as the main driver to read dump_atom file
 * ------------------------------------------------------------------ */
DumpAtom::DumpAtom(const char *infile, double *box, const double a_in, Crystal *crys) : Neigh(crys)
{
  FILE *fp = fopen(infile, "r");
  if (fp == NULL) return;

  a0 = a_in;

  for (int i=0; i<4; i++) fgets(oneline, MAXLINE, fp);
  natom = atoi(strtok(oneline," \t\n\r\f"));
  if (natom < 15) return;

  memory->create(x, natom+1,3,"DumpAtom:x");
  memory->create(attyp, natom+1,"DumpAtom:attyp");
  memory->create(neig_list, natom+1,NumNei,"DumpAtom:neig_list");

  fgets(oneline, MAXLINE, fp);
  double lo[3], hi[3], sh[3];
  sh[0] = sh[1] = sh[2] = 0.;
  for (int i=0; i<3; i++){
    fgets(oneline, MAXLINE, fp);
    lo[i] = atof(strtok(oneline," \t\n\r\f"));
    hi[i] = atof(strtok(NULL," \t\n\r\f"));
    char *ptr = strtok(NULL," \t\n\r\f");
    if (ptr) sh[i] = atof(ptr);
  }
  xlo = lo[0]; ylo = lo[1]; zlo = lo[2]; 
  xhi = hi[0]; yhi = hi[1]; zhi = hi[2]; 
  xy  = sh[0]; xz  = sh[1]; yz  = sh[2];
  box[0] = lx  = xhi-xlo;
  box[1] = ly  = yhi-ylo;
  box[2] = lz  = zhi-zlo;

  fgets(oneline, MAXLINE, fp);

  for (int i=0; i<natom; i++){
    fgets(oneline, MAXLINE, fp);
    int id = atoi(strtok(oneline," \t\n\r\f")); // get atom id
    int ip = atoi(strtok(NULL," \t\n\r\f"));
    attyp[id] = ip;
    x[id][0] = atof(strtok(NULL," \t\n\r\f"));  // get atomic position
    x[id][1] = atof(strtok(NULL," \t\n\r\f"));
    x[id][2] = atof(strtok(NULL," \t\n\r\f"));
  }
  fclose(fp);

  compute_neigh();

  flag = 0;

return;
}

/* ------------------------------------------------------------------
 * To supply the "local lattice constant" of atom id
 * ------------------------------------------------------------------ */
double DumpAtom::get_alat(int id)
{
  return a0;
}

/* ------------------------------------------------------------------
 * To free memory
 * ------------------------------------------------------------------ */
DumpAtom::~DumpAtom()
{
}

/* ------------------------------------------------------------------
 * Private method to compute the neighbor list for atoms with voro
 * neighbors less than NumNei, by searching all atoms
 * ------------------------------------------------------------------ */
void DumpAtom::compute_neigh()
{
  // assume BCC, get the nominal lattice constant
  double abcc = pow(lx*ly*lz/double(natom)*2., 1./3.);
  double cutoff2 = 3.0*abcc*abcc;

  int nmax = 50;
  int *count, **unsort_list;
  memory->create(count, natom+1, "compute_neigh:count");
  memory->create(unsort_list, nmax,natom+1,"compute_neigh:unsort_list");
  double **unsort_dist;
  memory->create(unsort_dist, nmax,natom+1,"compute:unsort_dist");
  for (int id=1; id<=natom; id++) count[id] = 0;

  // to get the neighbor list
  double Rij[3];
  for (int id=1; id<natom; id++){
    for (int jd=id+1; jd<=natom; jd++){
      for (int j=0; j<3; j++){
        Rij[j] = x[jd][j] - x[id][j];
        while (Rij[j] >  0.5) Rij[j] -= 1.;
        while (Rij[j] < -0.5) Rij[j] += 1.;
      }
      Rij[0] = Rij[0]*lx + Rij[1]*xy + Rij[2]*xz;
      Rij[1] = Rij[1]*ly + Rij[2]*yz;
      Rij[2] = Rij[2]*lz;
      double rij2 = Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2];

      if (rij2 >= cutoff2) continue;
      unsort_dist[count[id]][id] = -rij2; unsort_list[count[id]][id] = jd;
      unsort_dist[count[jd]][jd] = -rij2; unsort_list[count[jd]][jd] = id;

      count[id]++; count[jd]++;
      if (count[id] == nmax || count[jd] == nmax){
        nmax += 10;
        memory->grow(unsort_list,nmax,natom+1,"compute_neigh:unsort_list");
        memory->grow(unsort_dist,nmax,natom+1,"compute:unsort_dist");
      }
    }
  }

  for (int id=1; id<=natom; id++){
    local_prop.clear(); local_list.clear();
    for (int i=0; i<count[id]; i++){
      local_prop.push_back(unsort_dist[i][id]);
      local_list.push_back(unsort_list[i][id]);
    }
    int nlocal = local_list.size();
    if (nlocal < NumNei) printf("\nWarning: # of neighbors (%d) computed for %d is < NumNei!\n", nlocal, id);

    sort_list(id);

    // convert fractional coordinate into cartesian
    double xc[3];
    xc[0] = x[id][0] * lx + x[id][1] * xy + x[id][2] * xz;
    xc[1] = x[id][1] * ly + x[id][2] * yz;
    xc[2] = x[id][2] * lz;
    for (int j=0; j<3; j++) x[id][j] = xc[j];
  }
  memory->destroy(unsort_dist);
  memory->destroy(unsort_list);
  memory->destroy(count);

return;
}
