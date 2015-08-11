#ifndef CrysHCP_H
#define CrysHCP_H

#include "crystal.h"

class HCP : public Crystal {
public:

  HCP(const double);

  void align(double (*)[3]);
  double comp_alat(const double, const double);

private:
  double ca;
};
#endif
