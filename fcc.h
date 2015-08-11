#ifndef CrysFCC_H
#define CrysFCC_H

#include "crystal.h"

class FCC : public Crystal {
public:

  FCC();

  void align(double (*)[3]);
  double comp_alat(const double, const double);

};
#endif
