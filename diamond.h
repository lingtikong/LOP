#ifndef CrysDiamond_H
#define CrysDiamond_H

#include "crystal.h"

class Diamond : public Crystal {
public:

  Diamond();

  void align(double (*)[3]);
  double comp_alat(const double, const double);

};
#endif
