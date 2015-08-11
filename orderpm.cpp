#include "math.h"
#include "orderpm.h"
#include "timer.h"
#include "version.h"

/* ------------------------------------------------------------------
 * Constructor as the main driver to evaluate order paremter (rmsd)
 * ------------------------------------------------------------------ */
OrderPM::OrderPM(int narg, char **arg)
{
  voro = NULL;
  rmsd = NULL;
  crystal = NULL; fcrys_type = 0;
  fvoro = fout = frotated = NULL;

  opm = NULL;
  type = NULL;

  silent = 0;
  flag_type = 0;
  flag_rotated = 0;

  pbc[0] = pbc[1] = pbc[2] = 0;
  box[0] = box[1] = box[2] = 0.;

  int flag_comp_neig = 0;
  int ftype = 1;
  double a_in = 1.;
  double ca0 = sqrt(8./3.); // c/a ratio for hcp

  // analyze command line
  int iarg = 1;
  while (iarg < narg){
    if (strcmp(arg[iarg],"-o") == 0){  // define output file name
      if (++iarg >= narg) help();
      if (fout) delete []fout;
      fout = new char[strlen(arg[iarg])+1];
      strcpy(fout, arg[iarg]);

    } else if (strcmp(arg[iarg],"-p") == 0){ // define periodicity for PBC
      if (iarg+3 >= narg) help();
      box[0] = atof(arg[++iarg]);
      box[1] = atof(arg[++iarg]);
      box[2] = atof(arg[++iarg]);

    } else if (strcmp(arg[iarg], "-w") == 0){ // to output first reference local crystal and all rotated local crystals
      if (++iarg >= narg) help();
      if (frotated) delete []frotated;
      frotated = new char[strlen(arg[iarg])+1];
      strcpy(frotated, arg[iarg]);
      flag_rotated = 1;
      
    } else if (strcmp(arg[iarg], "-t") == 0){ // to define the threshold of order parameter distinguishing a crystal from a non-crystal
      if (++iarg >= narg) help();
      threshold = atof(arg[iarg]);
      if (threshold > 0.) flag_type = 1;

    } else if (strcmp(arg[iarg], "-i") == 0){ // to define the input file type: 1, voro++ custom output; 2, refined voro of dumpana; 3, lammps atom dump
      if (++iarg >= narg) help();       // by default: voro custom
      ftype = atoi(arg[iarg]);
      
    } else if (strcmp(arg[iarg], "-c") == 0){ // method to get neighbor list: 0, from voro; 1, compute if # from voro < NumNei; 2, compute;
      if (++iarg >= narg) help();       // in case the read file is atom style lammps dump, it carries the reference lattice constant.
      flag_comp_neig = atoi(arg[iarg]);
      a_in = atof(arg[iarg]);

    } else if (strcmp(arg[iarg], "-s") == 0){ // to work in silent mode, minimum output to standard output
      silent = 1;

    } else if (strcmp(arg[iarg], "-b") == 0){ // crystal will be bcc
      fcrys_type = 0;

    } else if (strcmp(arg[iarg], "-f") == 0){ // crystal will be fcc
      fcrys_type = 1;

    } else if (strcmp(arg[iarg], "-hcp") == 0){ // crystal will be hcp
      fcrys_type = 2;

    } else if (strcmp(arg[iarg], "-d") == 0){ // crystal will be hcp
      fcrys_type = 3;

    } else if (strcmp(arg[iarg], "-ca") == 0){ // to define the c/a ratio for hcp
      if (++iarg >= narg) help();
      ca0 = atof(arg[iarg]);

    } else if (strcmp(arg[iarg], "-h") == 0){ // to display help info
      help();

    } else { // mandatory input file; expected from voro++ or dump_atom from LAMMPS
      if (fvoro) delete []fvoro;
      fvoro = new char[strlen(arg[iarg])+1];
      strcpy(fvoro, arg[iarg]);
      
    }
    iarg++;
  }
  if (fvoro == NULL) help();

  // Show current code version info
  if (silent != 1) ShowVersion();
  // set default output file name if not supplied
  if (fout == NULL){
    fout = new char [9];
    strcpy(fout, "lopm.dat");
  }

  // initialize fractional reference local crystal and global methods
  // should be specified before read voro/dump info, as NumNei et al were set here
  if      (fcrys_type == 0) crystal = new BCC();
  else if (fcrys_type == 1) crystal = new FCC();
  else if (fcrys_type == 2) crystal = new HCP(ca0);
  else if (fcrys_type == 3) crystal = new Diamond();

  NumNei = crystal->NumNei;
  N1stN  = crystal->N1stN;
  memory = crystal->memory;

  // if voro++ custom, these are provided by command line option if required
  for (int idim=0; idim<3; idim++){
    if (box[idim] > 0.) pbc[idim] = 1;
    hbox[idim] = 0.5*box[idim];
  }

  Timer *timer;
  if (silent == 0) timer = new Timer();
  // get neighbor list and x info from input file
  if (ftype == 3){
    voro = new DumpAtom(fvoro, &box[0], a_in, crystal);

  } else if (ftype == 2){
    voro = new RefinedVoro(fvoro,&box[0],flag_comp_neig,crystal);

  } else {
    voro = new Voro(fvoro, box, flag_comp_neig, crystal);

  }

  for (int idim=0; idim<3; idim++){
    if (box[idim] > 0.) pbc[idim] = 1;
    hbox[idim] = 0.5*box[idim];
  }
  if (voro->flag){
    delete timer;
    printf("\nRead file %s failed. Program terminated.\n", fvoro);
    return;
  }

  // to compute the order parameter (rmsd)
  rmsd = new RMSD();
  memory->create(opm, voro->natom+1, "OrderPM:opm");
  compute_order_pm();

  // define atom type if a threshold is defined
  if (flag_type) compute_type();

  // to output the result
  output();

  if (silent == 0){
    timer->stop(); timer->print(); delete timer;
  }

return;
}

/* ------------------------------------------------------------------
 * Deconstructor to free memory
 * ------------------------------------------------------------------ */
OrderPM::~OrderPM()
{
  if (fvoro) delete []fvoro;
  if (fout)  delete []fout;
  if (frotated) delete []frotated;

  memory->destroy(opm);
  memory->destroy(type);

  memory = NULL;

  if (rmsd)  delete rmsd;
  if (voro)  delete voro;
  if (crystal) delete crystal;

return;
}

/* ------------------------------------------------------------------
 * Private method to compute the order parameter (rmsd)
 * ------------------------------------------------------------------ */
void OrderPM::compute_order_pm()
{

  if (flag_rotated){
    fpr = fopen(frotated, "w");
    write_1st_ref();
  }

  int natom = voro->natom;
  double pos_ref[NumNei+1][3], pos_loc[NumNei+1][3];

  for (int id=1; id<= natom; id++){

    double alat = voro->get_alat(id);
    // assign reference position
    for (int i=0; i<= NumNei; i++)
    for (int idim=0; idim<3; idim++) pos_ref[i][idim] = crystal->pos_ideal[i][idim]*alat;

    // assign neighbors of id to pos_loc
    for (int i=0; i<NumNei; i++){
      int jd = voro->neig_list[id][i];

      double Rij[3];
      for (int idim=0; idim<3; idim++){
        Rij[idim] = voro->x[jd][idim] - voro->x[id][idim];

        if (pbc[idim]){ // apply pbc
          while (Rij[idim] >= hbox[idim]) Rij[idim] -= box[idim];
          while (Rij[idim] < -hbox[idim]) Rij[idim] += box[idim];
        }

        pos_loc[i][idim] = Rij[idim];
      }
    }

    crystal->align(pos_loc);
    for (int idim=0; idim<3; idim++) pos_loc[NumNei][idim] = 0.;

    // get the order parameter by calling rmsd
    double mcom[3], v2ref[3];

    rmsd->calculate_rotation_rmsd(pos_ref, pos_loc, NumNei+1, mcom, v2ref, U, &opm[id]);

    if (flag_rotated) write_rotated(id, pos_loc);
  }
  if (flag_rotated) fclose(fpr);

return;
}

/* ------------------------------------------------------------------
 * Private method to identify the atom type based on the order
 * parameter computed and the threshold given.
 * ------------------------------------------------------------------ */
void OrderPM::compute_type()
{
  int natom = voro->natom;

  const int NMaxIt = 1000;
  int *swap, *otype;
  memory->create(otype, natom+1, "compute_type:otype");
  memory->create(type,  natom+1, "compute_type:type");

  for (int i=1; i<= natom; i++){
    if (opm[i] <= threshold) type[i] = 0; // crystal
    else type[i] = 2;                     // liquid or amorphous
  }

  // if nearly all neighbors of a crystal are liquid/amorphous, set it to be liquid/amorphous; and vice versa
  const int TwoNNei = NumNei+NumNei;
  int nit = 0;
  while (nit < NMaxIt){
    nit++;
    swap = otype; otype = type; type = swap; swap = NULL;

    int nreset = 0;
    for (int id=1; id<= natom; id++){
      int neisum = 0;
      for (int in=0; in<NumNei; in++){
        int jd = voro->neig_list[id][in];
        neisum += otype[jd];
      }

      type[id] = otype[id];
      if ( (otype[id] == 2) && (neisum <= 2        ) ){ type[id] = 0; nreset++; }
      if ( (otype[id] == 0) && (neisum >= TwoNNei-2) ){ type[id] = 2; nreset++; }
    }
    if (nreset == 0) break;
  }
  if (nit >= NMaxIt && silent==0 ) printf("\nWarning: your cutoff might be improper!\n\n");

  swap = otype; otype = type; type = swap; swap = NULL;
  // assign type: 0, crystal; 1, crystal at interface; 2, liquid/amorphous; 3, liquid/amorhpous at interface
  neach[0] = neach[1] = neach[2] = neach[3] = 0;
  for (int id=1; id<=natom; id++){
    int it = otype[id];
    int neisum = 0;
    for (int in=0; in<NumNei; in++){
      int jd = voro->neig_list[id][in];
      neisum += otype[jd];
    }
    if (neisum > 4 && neisum<TwoNNei-4) it++;
    else it = (neisum+4)/NumNei;

    type[id] = it;
  }

  nit = 0;
  while (nit < NMaxIt){
    nit++;
    swap = otype; otype = type; type = swap; swap = NULL;
    neach[0] = neach[1] = neach[2] = neach[3] = 0;

    int nreset = 0;
    for (int id=1; id<=natom; id++){
      int has[4];
      has[0] = has[1] = has[2] = has[3] = 0;
      for (int in=0; in<NumNei; in++){
        int jd = voro->neig_list[id][in];
        int jt = otype[jd];
        has[jt] = 1;
      }

      int it = otype[id];
      type[id] = it;

      if (it == 0 && has[2] == 1) { type[id] = 1; nreset++; }
      if (it == 2 && has[0] == 1) { type[id] = 3; nreset++; }
      if (it == 1){
        if(has[0] == 0){
          nreset++;
          if (has[3]) type[id] = 3;
          else type[id] = 2;
        }
        if (has[2] == 0 && has[3] == 0){
          nreset++;
          type[id] = 0;
        }
      }
      if (it == 3){
        if (has[2] == 0){
          nreset++;
          if (has[1]) type[id] = 1;
          else type[id] = 0;
        }
        if (has[0] == 0 && has[1] == 0){
          nreset++;
          type[id] = 2;
        }
      }

      neach[type[id]]++;
    }

    if (nreset == 0) break;
  }
  if (nit >= NMaxIt && silent==0) printf("\nWarning: your cutoff might be inappropriate!\n\n");

  memory->destroy(otype);

return;
}

/* ------------------------------------------------------------------
 * Private method to write out the result
 * ------------------------------------------------------------------ */
void OrderPM::output()
{
  int natom = voro->natom;
  double **x = voro->x;

  FILE *fp = fopen(fout, "w");
  fprintf(fp,"# box info: %lg %lg %lg\n", box[0], box[1], box[2]);
  if (flag_type){
    fprintf(fp,"# threshold for type identification: %lg\n", threshold);
    fprintf(fp,"# Type: 0, crystal;     1, crystal at interface;\n#\t2, non-crystal; 3, interfacial non-crystal.\n");
    fprintf(fp,"# Num for each type: %d %d %d %d\n#\n", neach[0], neach[1], neach[2], neach[3]);
    if (voro->attyp){
      fprintf(fp,"# id\tx\ty\tz\topm type\tattyp\n");
      for (int id=1; id<= natom; id++)
      fprintf(fp,"%d %lg %lg %lg %lg %d %d\n", id, x[id][0], x[id][1], x[id][2], opm[id], type[id],voro->attyp[id]);
    } else {
      fprintf(fp,"# id\tx\ty\tz\topm type\n");
      for (int id=1; id<= natom; id++)
      fprintf(fp,"%d %lg %lg %lg %lg %d\n", id, x[id][0], x[id][1], x[id][2], opm[id], type[id]);
    }

  } else {
    if (voro->attyp){
      fprintf(fp,"# id\tx\ty\tz\topm\tattyp\n");
      for (int id=1; id<= natom; id++)
      fprintf(fp,"%d %lg %lg %lg %lg %d\n", id, x[id][0], x[id][1], x[id][2], opm[id],voro->attyp[id]);

    } else {
      fprintf(fp,"# id\tx\ty\tz\topm\n");
      for (int id=1; id<= natom; id++)
      fprintf(fp,"%d %lg %lg %lg %lg\n", id, x[id][0], x[id][1], x[id][2], opm[id]);
    }
  }
  fclose(fp);

return;
}

/* ------------------------------------------------------------------
 * Private method to write one rotated local crystal
 * ------------------------------------------------------------------ */
void OrderPM::write_rotated(int id, double lpos[][3])
{
  double pos_rotated[NumNei+1][3];
  for (int i=0; i<= NumNei; i++){
    for (int j=0; j<3; j++){
      pos_rotated[i][j] = 0.;
      for (int k=0; k<3; k++) pos_rotated[i][j] += U[j][k]*lpos[i][k];
    }
  }

  int ip = 1;
  fprintf(fpr,"%d\n", NumNei+1);
  fprintf(fpr,"Local cluster center on atom id = %d, rmsd= %lg\n", id, opm[id]);
  for (int i=0; i< NumNei; i++){
    int jd = voro->neig_list[id][i];
    if (voro->attyp) ip = voro->attyp[jd];
    fprintf(fpr,"%d %lg %lg %lg %d\n", ip, pos_rotated[i][0], pos_rotated[i][1], pos_rotated[i][2], jd);
  }
  if (voro->attyp) ip = voro->attyp[id];
  fprintf(fpr,"%d %lg %lg %lg %d\n", ip, pos_rotated[NumNei][0], pos_rotated[NumNei][1], pos_rotated[NumNei][2], id);

return;
}

/* ------------------------------------------------------------------
 * Private method to write the reference local crystal; the lattice
 * constant is given by the "local lattice constant" of the first atom
 * ------------------------------------------------------------------ */
void OrderPM::write_1st_ref()
{
  double pos_ref[NumNei+1][3];
  // assign reference position
  for (int i=0; i<= NumNei; i++)
  for (int j=0; j<3; j++) pos_ref[i][j] = crystal->pos_ideal[i][j]*voro->get_alat(1);

  fprintf(fpr,"%d\n", NumNei+1);
  fprintf(fpr,"Reference cluster for the first atom.\n");
  for (int i=0; i<= NumNei; i++) fprintf(fpr,"%d %lg %lg %lg\n", i/NumNei+1, pos_ref[i][0], pos_ref[i][1], pos_ref[i][2]);
}

/* ------------------------------------------------------------------
 * Private method to print out help info
 * ------------------------------------------------------------------ */
void OrderPM::help()
{
  ShowVersion();
  for (int i=0; i<21 ;i++) printf("----");
  printf("\nUsage: lop [options] file\n\n");
  printf("  file is either:\n");
  printf("    1) output from voro++ with custom output options: %%i %%q %%v %%s %%F %%f %%n\".\n");
  printf("    2) output from refined voro computation of dumpana;\n");
  printf("    3) lammps atom style dump file;\n");
  printf("Available options:\n");
  printf("  -o outfile       To define the output file name; default: lopm.dat\n");
  printf("  -p lx ly lz      To define the box size for PBC, a non-positive number suggest no pbc; default: 0\n");
  printf("                   Only necessary if file is read from voro++ custom output.\n");
  printf("  -w rotated_file  Ask for output the rotated local clusters to file; default: no output.\n");
  printf("  -t threshold     To define the threshold to separate a crystal from a liquid/amorphous. default: not set.\n");
  printf("  -i filetype      To define the input file type: 1, voro++ custom output; 2, voro from dumpana; 3, lammps atom dump.\n");
  printf("                     Default: 1\n");
  printf("  -c method        To define the way to get neighbor list; 0, from voro file; 1, compute only if # of neighbors\n");
  printf("                     from voro info is less than NumNei; 2, compute from atomic positions read from voro.\n");
  printf("                     In case file is lammps atom dump, the value of method is acutally the reference lattice constant.\n");
  printf("                     Default: 0 for filetyle != 3, 1 for filetype == 3.\n");
  printf("  -s               To work in silent mode, no warning or timing output. Default: not set.\n");
  printf("  -b               To define the reference lattice as bcc; default.\n");
  printf("  -f               To define the reference lattice as fcc;\n");
  printf("  -hcp             To define the reference lattice as hcp;\n");
  printf("  -d               To define the reference lattice as diamond (only 4 neighbors);\n");
  printf("  -ca c/a          To define the c/a ratio for hcp, default:sqrt(8/3).\n");
  printf("  -h               To display this help info.\n\n");
  exit(0);
return;
}

/* ------------------------------------------------------------------
 * Private method to show code version
 * ------------------------------------------------------------------ */
void OrderPM::ShowVersion()
{
  printf("\n                _      ____  _____  \n");
  printf("               | |    / __ \\|  __ \\ \n");
  printf("               | |   | |  | | |__) |  \n");
  printf("               | |   | |  | |  ___/   \n");
  printf("               | |___| |__| | |       \n");
  printf("               |______\\____/|_|      \n");
  printf("\nCode to compute the Local Order Parameter, defined as the RMSD with respect to\n");
  printf("the respective ideal crystal. Version 1.%d, compiled on %s.\n", VERSION, __DATE__);

return;
}

/* ------------------------------------------------------------------ */
