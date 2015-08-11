#ifndef VORO_H
#define VORO_H

#include "neigh.h"

class Voro : public Neigh {
public:
  Voro(const char *,const double [],const int, Crystal *);
  ~Voro();

  double get_alat(int);

private:
  int pbc[3];
  double prd[3], hprd[3];
  int flag_comp_neigh;
  void compute_neigh();
  std::vector<int> atom_list;

  std::vector<double> a0, vol;

};
#endif
