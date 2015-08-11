#ifndef REFINED_VORO_H
#define REFINED_VORO_H

#include "neigh.h"

class RefinedVoro : public Neigh {
public:
  RefinedVoro(const char *, double *,const int, Crystal *);
  ~RefinedVoro();

  double get_alat(int);

private:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double lx, ly, lz;
  double xy, xz, yz;
  double hbox[3], fbox[3];

  int flag_comp_neigh;
  void compute_neigh();
  std::vector<int> atom_list;

  std::vector<double> a0, vol;

};
#endif
