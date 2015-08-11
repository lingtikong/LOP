#ifndef CrysBCC_H
#define CrysBCC_H

#include "crystal.h"

class BCC : public Crystal {
public:

  BCC();

  void align(double (*)[3]);
  double comp_alat(const double, const double);

};
#endif
