#ifndef Crystal_H
#define Crystal_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "math.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

class Crystal {
public:

  Crystal();
  ~Crystal();

  virtual void align(double (*)[3]){return;};
  virtual double comp_alat(const double, const double){return 1.;};

  int NumNei, N1stN;
  double **pos_ideal;

  Memory *memory;

  double angle(double *, double *);
  double norm2(double *);
  void cross(double *, double *, double *); // 1 = 2 x 3
  void normalize(double *);

};
#endif
