#ifndef DUMP_ATOM_H
#define DUMP_ATOM_H

#include "neigh.h"

class DumpAtom : public Neigh {
public:
  DumpAtom(const char *, double *, const double, Crystal *);
  ~DumpAtom();

  double get_alat(int);

private:
  double a0;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double lx, ly, lz;
  double xy, xz, yz;

  void compute_neigh();

};
#endif
