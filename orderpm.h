#ifndef ORDER_PARAMETER_H
#define ORDER_PARAMETER_H

#include "neigh.h"
#include "voro.h"
#include "dump_atom.h"
#include "refinedvoro.h"
#include "rmsd.h"
#include "lists.h"

using namespace std;

class OrderPM {
public:
  OrderPM(int, char **);
  ~OrderPM();

private:
  Neigh *voro;                 // neighbor info: x, neigh, alat
  RMSD *rmsd;                  // method to get rmsd

  Crystal *crystal;
  Memory *memory;
  int fcrys_type;              // flag to indicate crystal type: 0, bcc; 1, fcc
  int NumNei, N1stN;

  int pbc[3];
  char *fvoro, *fout;
  double box[3], hbox[3];

  int flag_type, silent;
  double threshold;
  double *opm;                 // computed order paremter
  int *type, neach[4];         // 0: crystal; 1: interfacial crystal; 2: non-crystal; 3: interfacial non-crystal

  double U[3][3];

  FILE *fpr;
  char *frotated;
  int flag_rotated;
  void write_rotated(int, double (*)[3]);
  void write_1st_ref();

  void compute_order_pm();
  void compute_type();
  void output();

  void help();
  void ShowVersion();
};

#endif
